/**
 * UGENE - Integrated Bioinformatics Tools.
 * Copyright (C) 2008-2026 UniPro <ugene@unipro.ru>
 * http://ugene.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */

#include "RestrictionMapWidget.h"

#include <climits>
#include <atomic>
#include <memory>
#include <utility>

#include <QButtonGroup>
#include <QElapsedTimer>
#include <QHash>
#include <QHeaderView>
#include <QHBoxLayout>
#include <QMutex>
#include <QPainter>
#include <QPixmap>
#include <QPushButton>
#include <QRunnable>
#include <QSet>
#include <QSizePolicy>
#include <QThreadPool>
#include <QTimer>
#include <QVBoxLayout>

#include <U2Algorithm/EnzymeModel.h>

#include <U2Core/Annotation.h>
#include <U2Core/AnnotationGroup.h>
#include <U2Core/AnnotationSelection.h>
#include <U2Core/AnnotationTableObject.h>
#include <U2Core/AppContext.h>
#include <U2Core/AutoAnnotationsSupport.h>
#include <U2Core/DNAAlphabet.h>
#include <U2Core/DNASequenceObject.h>
#include <U2Core/DNASequenceSelection.h>
#include <U2Core/Settings.h>
#include <U2Core/U1AnnotationUtils.h>
#include <U2Core/U2OpStatusUtils.h>

#include <U2Gui/GUIUtils.h>

#include <U2View/ADVSequenceObjectContext.h>
#include <U2View/AnnotatedDNAView.h>
#include <U2View/AutoAnnotationUtils.h>

#include "HostMethylation.h"

#define ENZYME_FOLDER_ITEM_TYPE 1022
#define ENZYME_ITEM_TYPE 1023

namespace {

const char* SHOW_EMPTY_ENZYMES_KEY = "circular_view/restriction_map_show_empty_enzymes";
const int UNLIMITED_HIT_COUNT = INT_MAX;
// Lists this size stay on the calling turn. Larger updates yield so the cursor can move.
// Enzyme rows past this size keep a count and build child rows when the row is expanded.
constexpr int INTERACTIVE_TREE_LIMIT = 200;
constexpr int INTERACTIVE_SLICE_MS = 12;
constexpr int METHYLATION_COLUMN_WIDTH = 180;
constexpr int MATERIALIZED_SITES_PER_ENZYME = 250;

QIcon blockedMethylationIcon() {
    static QIcon icon;
    if (!icon.isNull()) {
        return icon;
    }
    QPixmap pixmap(32, 32);
    pixmap.setDevicePixelRatio(2.0);
    pixmap.fill(Qt::transparent);
    QPainter painter(&pixmap);
    painter.setRenderHint(QPainter::Antialiasing, true);
    QPen pen(QColor(0xC4, 0x16, 0x16), 2.2, Qt::SolidLine, Qt::RoundCap);
    painter.setPen(pen);
    painter.drawLine(QPointF(3.5, 3.5), QPointF(12.5, 12.5));
    painter.drawLine(QPointF(12.5, 3.5), QPointF(3.5, 12.5));
    painter.end();
    icon = QIcon(pixmap);
    return icon;
}

QString methylationLabel(int flags, int damCount, int dcmCount) {
    QStringList words;
    if (flags & U2::HostMethylationDam) {
        words << (damCount < 0 ? QStringLiteral("Dam") : QStringLiteral("Dam (%1)").arg(damCount));
    }
    if (flags & U2::HostMethylationDcm) {
        words << (dcmCount < 0 ? QStringLiteral("Dcm") : QStringLiteral("Dcm (%1)").arg(dcmCount));
    }
    return words.join(QStringLiteral(" "));
}

QString blockedSiteCountToolTip(int count, bool dam) {
    if (count == 1) {
        return dam ? U2::RestrctionMapWidget::tr("1 site is blocked by Dam methylation in a Dam+ host")
                   : U2::RestrctionMapWidget::tr("1 site is blocked by Dcm methylation in a Dcm+ host");
    }
    return dam ? U2::RestrctionMapWidget::tr("%1 sites are blocked by Dam methylation in a Dam+ host").arg(count)
               : U2::RestrctionMapWidget::tr("%1 sites are blocked by Dcm methylation in a Dcm+ host").arg(count);
}

QString folderMethylationToolTip(int damCount, int dcmCount) {
    QStringList parts;
    if (damCount > 0) {
        parts << blockedSiteCountToolTip(damCount, true);
    }
    if (dcmCount > 0) {
        parts << blockedSiteCountToolTip(dcmCount, false);
    }
    return parts.join(QStringLiteral("\n"));
}

void applyMethylationMarker(QTreeWidgetItem* item, int flags, int damCount, int dcmCount) {
    if (flags == U2::HostMethylationNone) {
        item->setIcon(1, QIcon());
        item->setText(1, QString());
        item->setToolTip(1, QString());
        item->setData(1, Qt::ForegroundRole, QVariant());
        return;
    }
    item->setIcon(1, blockedMethylationIcon());
    item->setText(1, methylationLabel(flags, damCount, dcmCount));
    if (damCount < 0 && dcmCount < 0) {
        QStringList parts;
        if (flags & U2::HostMethylationDam) {
            parts << U2::RestrctionMapWidget::tr("Blocked by Dam methylation in a Dam+ host");
        }
        if (flags & U2::HostMethylationDcm) {
            parts << U2::RestrctionMapWidget::tr("Blocked by Dcm methylation in a Dcm+ host");
        }
        item->setToolTip(1, parts.join(QStringLiteral("\n")));
    } else {
        item->setToolTip(1, folderMethylationToolTip(damCount, dcmCount));
    }
    item->setForeground(1, QColor(0xC4, 0x16, 0x16));
}

}  // namespace

namespace U2 {

//////////////////////////////////////////////////////////////////////////
/// EnzymeItem

EnzymeItem::EnzymeItem(const QString& location, Annotation* a, quint64 itemToken)
    : QTreeWidgetItem(QStringList(location), ENZYME_ITEM_TYPE), annotation(a), token(itemToken) {
}

void EnzymeItem::setMethylationFlags(int flags) {
    methylationFlags = flags;
    applyMethylationMarker(this, flags, -1, -1);
}

//////////////////////////////////////////////////////////////////////////
/// EnzymeFolderItem

EnzymeFolderItem::EnzymeFolderItem(const QString& name)
    : QTreeWidgetItem(ENZYME_FOLDER_ITEM_TYPE), enzymeName(name) {
    setText(0, QString("%1 : %2 %3").arg(name).arg(0).arg("sites"));
}

EnzymeItem* EnzymeFolderItem::createEnzymeItem(Annotation* enzAnn, quint64 token) {
    const SharedAnnotationData& data = enzAnn->getData();
    const QString location = U1AnnotationUtils::buildLocationString(data);
    return new EnzymeItem(location, enzAnn, token);
}

void EnzymeFolderItem::detachEnzymeItem(EnzymeItem* item) {
    adjustBlockedCounts(item->getMethylationFlags(), HostMethylationNone);
    removeChild(item);
    delete item;
}

void EnzymeFolderItem::publishSiteCount() {
    const int count = siteCount();
    setIcon(0, QIcon(count == 0 ? ":circular_view/images/empty_folder.png" : ":circular_view/images/folder.png"));
    const QString site = count == 1 ? RestrctionMapWidget::tr("site") : RestrctionMapWidget::tr("sites");
    setText(0, QString("%1 : %2 %3").arg(getName()).arg(count).arg(site));
}

void EnzymeFolderItem::addHeldSite() {
    ++heldSiteCount;
    syncExpander();
}

void EnzymeFolderItem::removeHeldSite() {
    if (heldSiteCount > 0) {
        --heldSiteCount;
    }
    syncExpander();
}

void EnzymeFolderItem::syncExpander() {
    setChildIndicatorPolicy(heldSiteCount > 0 ? QTreeWidgetItem::ShowIndicator : QTreeWidgetItem::DontShowIndicatorWhenChildless);
}

void EnzymeFolderItem::adjustBlockedCounts(int oldFlags, int newFlags) {
    if ((oldFlags & HostMethylationDam) != 0 && (newFlags & HostMethylationDam) == 0) {
        --damBlockedCount;
    } else if ((oldFlags & HostMethylationDam) == 0 && (newFlags & HostMethylationDam) != 0) {
        ++damBlockedCount;
    }
    if ((oldFlags & HostMethylationDcm) != 0 && (newFlags & HostMethylationDcm) == 0) {
        --dcmBlockedCount;
    } else if ((oldFlags & HostMethylationDcm) == 0 && (newFlags & HostMethylationDcm) != 0) {
        ++dcmBlockedCount;
    }
}

void EnzymeFolderItem::publishMethylationMarker() {
    int flags = HostMethylationNone;
    if (damBlockedCount > 0) {
        flags |= HostMethylationDam;
    }
    if (dcmBlockedCount > 0) {
        flags |= HostMethylationDcm;
    }
    applyMethylationMarker(this, flags, damBlockedCount, dcmBlockedCount);
}

//////////////////////////////////////////////////////////////////////////
/// RestrictionMapWidget

struct RestrctionMapWidget::MethylationState {
    struct Site {
        quint64 token = 0;
        QString enzymeName;
        QVector<U2Region> regions;
    };
    struct Mark {
        quint64 token = 0;
        int flags = 0;
    };
    struct HeldSite {
        Annotation* annotation = nullptr;
        int flags = 0;
    };
    // The worker publishes marks here. The widget drops its pointer to cancel.
    struct Mailbox {
        std::atomic<bool> cancelled{false};
        QMutex mutex;
        bool ready = false;
        QVector<Mark> marks;
    };

    int epoch = 0;
    quint64 nextToken = 1;
    bool insertRunning = false;
    bool insertScheduled = false;
    bool detachRunning = false;
    bool detachScheduled = false;
    bool jobActive = false;
    bool pollScheduled = false;
    bool applyScheduled = false;
    bool columnLocked = false;
    bool bulkFill = false;
    bool materializeScheduled = false;
    bool sequenceReady = false;
    bool circular = false;
    QByteArray sequence;

    QVector<Annotation*> pendingAnnotations;
    int pendingAnnotationCursor = 0;
    QSet<Annotation*> pendingAnnotationLive;
    QHash<Annotation*, EnzymeItem*> itemsByAnnotation;
    QHash<quint64, EnzymeItem*> itemsByToken;
    QHash<QString, EnzymeFolderItem*> foldersByName;
    QHash<quint64, HeldSite> heldByToken;
    QHash<Annotation*, quint64> heldTokenByAnnotation;
    QHash<QString, QVector<quint64>> heldTokensByEnzyme;
    QVector<quint64> pendingMaterialize;
    int pendingMaterializeCursor = 0;
    QVector<EnzymeItem*> pendingDetach;
    QVector<Site> pendingSites;
    QVector<Mark> pendingMarks;
    int pendingMarkCursor = 0;
    std::shared_ptr<Mailbox> mailbox;
};

RestrctionMapWidget::RestrctionMapWidget(ADVSequenceObjectContext* context, QWidget* p)
    : QWidget(p), ctx(nullptr), annotatedDnaView(nullptr), treeWidget(nullptr), methylation(std::make_unique<MethylationState>()) {
    assert(context != nullptr);
    setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Expanding);
    annotatedDnaView = context->getAnnotatedDNAView();

    auto layout = new QVBoxLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);

    treeWidget = new QTreeWidget(this);
    treeWidget->setObjectName("restrictionMapTreeWidget");
    treeWidget->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Expanding);
    treeWidget->setColumnCount(2);
    treeWidget->setHeaderLabels({tr("Restriction Sites Map"), QString()});
    treeWidget->header()->setStretchLastSection(false);
    treeWidget->header()->setSectionResizeMode(0, QHeaderView::Stretch);
    treeWidget->header()->setSectionResizeMode(1, QHeaderView::ResizeToContents);
    connect(treeWidget, SIGNAL(itemSelectionChanged()), SLOT(sl_itemSelectionChanged()));
    connect(treeWidget, SIGNAL(itemExpanded(QTreeWidgetItem*)), SLOT(sl_itemExpanded(QTreeWidgetItem*)));

    layout->addWidget(treeWidget, 1);
    initFooterButtons();

    connect(annotatedDnaView, SIGNAL(si_activeSequenceWidgetChanged(ADVSequenceWidget*, ADVSequenceWidget*)), SLOT(sl_onActiveSequenceChanged()));
    connect(annotatedDnaView, SIGNAL(si_sequenceRemoved(ADVSequenceObjectContext*)), SLOT(sl_onSequenceRemoved(ADVSequenceObjectContext*)));
    setSequenceContext(context);
}

RestrctionMapWidget::~RestrctionMapWidget() {
    invalidateMethylationTracking();
}

void RestrctionMapWidget::setSequenceContext(ADVSequenceObjectContext* sequenceContext) {
    if (ctx == sequenceContext) {
        return;
    }
    unregisterAnnotationObjects();
    ctx = sequenceContext;
    if (ctx == nullptr) {
        invalidateMethylationTracking();
        treeWidget->clear();
        return;
    }
    connect(ctx, SIGNAL(si_annotationObjectAdded(AnnotationTableObject*)), SLOT(sl_onAnnotationObjectAdded(AnnotationTableObject*)));
    connect(ctx, SIGNAL(si_annotationObjectRemoved(AnnotationTableObject*)), SLOT(sl_onAnnotationObjectRemoved(AnnotationTableObject*)));
    registerAnnotationObjects();
    rebuildTree();
}

void RestrctionMapWidget::unregisterAnnotationObjects() {
    if (ctx == nullptr) {
        return;
    }
    QSet<AnnotationTableObject*> annotationObjects = ctx->getAnnotationObjects(true);
    foreach (AnnotationTableObject* annotationObject, annotationObjects) {
        annotationObject->disconnect(this);
    }
    ctx->disconnect(this);
}

void RestrctionMapWidget::connectAnnotationObject(AnnotationTableObject* object) {
    connect(object, SIGNAL(si_onAnnotationsAdded(const QList<Annotation*>&)), SLOT(sl_onAnnotationsAdded(const QList<Annotation*>&)));
    connect(object, SIGNAL(si_onAnnotationsRemoved(const QList<Annotation*>&)), SLOT(sl_onAnnotationsRemoved(const QList<Annotation*>&)));
    connect(object, SIGNAL(si_onAnnotationsInGroupRemoved(const QList<Annotation*>&, AnnotationGroup*)), SLOT(sl_onAnnotationsInGroupRemoved(const QList<Annotation*>&, AnnotationGroup*)));
    connect(object, SIGNAL(si_onGroupCreated(AnnotationGroup*)), SLOT(sl_onAnnotationsGroupCreated(AnnotationGroup*)));
}

void RestrctionMapWidget::rebuildTree() {
    treeWidget->blockSignals(true);
    updateTreeWidget();
    initTreeWidget();
    treeWidget->blockSignals(false);
    updateEmptyEnzymeVisibility();
}

void RestrctionMapWidget::initFooterButtons() {
    // GroupOptionsWidget moves this child onto the row beside the Help button.
    auto footer = new QWidget(this);
    footer->setObjectName("optionsPanelFooter");
    auto footerLayout = new QHBoxLayout(footer);
    footerLayout->setContentsMargins(0, 0, 0, 0);
    footerLayout->setSpacing(4);

    auto hitGroup = new QButtonGroup(footer);
    hitGroup->setExclusive(true);
    const auto addHitButton = [&](const QString& text, const QString& objectName, const QString& toolTip) {
        auto button = new QPushButton(text, footer);
        button->setObjectName(objectName);
        button->setToolTip(toolTip);
        button->setCheckable(true);
        button->setFocusPolicy(Qt::NoFocus);
        button->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        hitGroup->addButton(button);
        footerLayout->addWidget(button);
        connect(button, &QPushButton::clicked, this, &RestrctionMapWidget::sl_maxHitsClicked);
        return button;
    };
    maxHits1Button = addHitButton(tr("1"), "restrictionMapHits1Button", tr("Find restriction sites with a maximum of 1 hit"));
    maxHits2Button = addHitButton(tr("2"), "restrictionMapHits2Button", tr("Find restriction sites with a maximum of 2 hits"));
    maxHitsAllButton = addHitButton(tr("All"), "restrictionMapHitsAllButton", tr("Find restriction sites with no maximum hit count"));

    Settings* settings = AppContext::getSettings();
    const bool limitHits = settings->getValue(EnzymeSettings::ENABLE_HIT_COUNT, false).toBool();
    const int maxHits = settings->getValue(EnzymeSettings::MAX_HIT_VALUE, UNLIMITED_HIT_COUNT).toInt();
    if (limitHits && maxHits == 1) {
        maxHits1Button->setChecked(true);
    } else if (limitHits && maxHits == 2) {
        maxHits2Button->setChecked(true);
    } else if (!limitHits || maxHits >= UNLIMITED_HIT_COUNT) {
        maxHitsAllButton->setChecked(true);
    }

    footerLayout->addSpacing(8);
    auto showEmptyButton = new QPushButton(tr("Empty"), footer);
    showEmptyButton->setObjectName("restrictionMapShowEmptyButton");
    showEmptyButton->setToolTip(tr("Show enzymes with no sites. They stay in the list and are grayed out."));
    showEmptyButton->setCheckable(true);
    showEmptyButton->setFocusPolicy(Qt::NoFocus);
    showEmptyEnzymes = settings->getValue(SHOW_EMPTY_ENZYMES_KEY, true).toBool();
    showEmptyButton->setChecked(showEmptyEnzymes);
    connect(showEmptyButton, &QPushButton::toggled, this, &RestrctionMapWidget::sl_showEmptyEnzymesToggled);
    footerLayout->addWidget(showEmptyButton);
}

int RestrctionMapWidget::selectedMaxHitCount() const {
    if (maxHits1Button != nullptr && maxHits1Button->isChecked()) {
        return 1;
    }
    if (maxHits2Button != nullptr && maxHits2Button->isChecked()) {
        return 2;
    }
    return UNLIMITED_HIT_COUNT;
}

void RestrctionMapWidget::sl_maxHitsClicked() {
    const int maxHits = selectedMaxHitCount();
    const bool limitHits = maxHits != UNLIMITED_HIT_COUNT;
    Settings* settings = AppContext::getSettings();
    settings->setValue(EnzymeSettings::ENABLE_HIT_COUNT, limitHits);
    settings->setValue(EnzymeSettings::MIN_HIT_VALUE, 1);
    settings->setValue(EnzymeSettings::MAX_HIT_VALUE, maxHits);
    if (ctx != nullptr) {
        AutoAnnotationUtils::triggerAutoAnnotationsUpdate(ctx, ANNOTATION_GROUP_ENZYME);
    }
}

void RestrctionMapWidget::sl_showEmptyEnzymesToggled(bool showEmpty) {
    showEmptyEnzymes = showEmpty;
    AppContext::getSettings()->setValue(SHOW_EMPTY_ENZYMES_KEY, showEmpty);
    updateEmptyEnzymeVisibility();
}

void RestrctionMapWidget::updateEmptyEnzymeVisibility() {
    const QColor emptyText = treeWidget->palette().color(QPalette::Disabled, QPalette::Text);
    const int count = treeWidget->topLevelItemCount();
    for (int i = 0; i < count; ++i) {
        auto item = static_cast<EnzymeFolderItem*>(treeWidget->topLevelItem(i));
        const bool empty = item->siteCount() == 0;
        const bool hide = empty && !showEmptyEnzymes;
        item->setHidden(hide);
        if (hide && item->isSelected()) {
            item->setSelected(false);
        }
        // An explicit color sticks through selection, so only empty rows get the disabled gray.
        if (empty) {
            item->setForeground(0, emptyText);
        } else {
            item->setData(0, Qt::ForegroundRole, QVariant());
        }
    }
}

void RestrctionMapWidget::sl_onActiveSequenceChanged() {
    ADVSequenceObjectContext* activeContext = annotatedDnaView->getActiveSequenceContext();
    if (activeContext == nullptr || activeContext == ctx) {
        return;
    }
    const DNAAlphabet* alphabet = activeContext->getAlphabet();
    if (alphabet == nullptr || !alphabet->isNucleic()) {
        return;
    }
    setSequenceContext(activeContext);
}

void RestrctionMapWidget::sl_onSequenceRemoved(ADVSequenceObjectContext* sequenceContext) {
    if (sequenceContext != ctx) {
        return;
    }
    ADVSequenceObjectContext* replacement = nullptr;
    const QList<ADVSequenceObjectContext*> contexts = annotatedDnaView->getSequenceContexts();
    foreach (ADVSequenceObjectContext* candidate, contexts) {
        if (candidate == sequenceContext) {
            continue;
        }
        const DNAAlphabet* alphabet = candidate->getAlphabet();
        if (alphabet != nullptr && alphabet->isNucleic()) {
            replacement = candidate;
            break;
        }
    }
    setSequenceContext(replacement);
}

void RestrctionMapWidget::sl_onAnnotationObjectAdded(AnnotationTableObject* object) {
    connectAnnotationObject(object);
    sl_onAnnotationsAdded(object->getAnnotations());
}

void RestrctionMapWidget::sl_onAnnotationObjectRemoved(AnnotationTableObject* object) {
    object->disconnect(this);
    if (ctx != nullptr) {
        rebuildTree();
    }
}

void RestrctionMapWidget::updateTreeWidget() {
    invalidateMethylationTracking();
    treeWidget->clear();

    QString selection = AppContext::getSettings()->getValue(EnzymeSettings::LAST_SELECTION).toString();
    if (selection.isEmpty()) {
        selection = EnzymeSettings::COMMON_ENZYMES;
    }
    QStringList selectedEnzymes = selection.split(ENZYME_LIST_SEPARATOR, Qt::SkipEmptyParts);

    QList<QTreeWidgetItem*> items;
    foreach (const QString& enzyme, selectedEnzymes) {
        auto item = new EnzymeFolderItem(enzyme);
        item->setIcon(0, QIcon(":circular_view/images/empty_folder.png"));
        if (!methylation->foldersByName.contains(enzyme)) {
            methylation->foldersByName.insert(enzyme, item);
        }
        items.append(item);
    }
    treeWidget->insertTopLevelItems(0, items);
    treeWidget->sortItems(0, Qt::AscendingOrder);
    updateEmptyEnzymeVisibility();
}

void RestrctionMapWidget::registerAnnotationObjects() {
    QSet<AnnotationTableObject*> aObjs = ctx->getAnnotationObjects(true);
    foreach (AnnotationTableObject* ao, aObjs) {
        connectAnnotationObject(ao);
    }
}

bool RestrctionMapWidget::loadSequence(QByteArray& sequence, bool& circular) const {
    sequence.clear();
    circular = false;
    if (ctx == nullptr || ctx->getSequenceObject() == nullptr) {
        return false;
    }
    U2OpStatusImpl status;
    sequence = ctx->getSequenceData(U2Region(0, ctx->getSequenceLength()), status);
    if (status.hasError()) {
        sequence.clear();
        return false;
    }
    circular = ctx->getSequenceObject()->isCircular();
    return true;
}

void RestrctionMapWidget::sl_onAnnotationsAdded(const QList<Annotation*>& anns) {
    queueAnnotations(anns);
}

void RestrctionMapWidget::sl_onAnnotationsRemoved(const QList<Annotation*>& anns) {
    removeAnnotationsFromMap(anns);
}

EnzymeFolderItem* RestrctionMapWidget::findEnzymeFolderByName(const QString& enzymeName) {
    if (methylation == nullptr) {
        return nullptr;
    }
    return methylation->foldersByName.value(enzymeName, nullptr);
}

void RestrctionMapWidget::sl_itemSelectionChanged() {
    if (ctx == nullptr) {
        return;
    }
    Annotation* selectedAnnotation = nullptr;
    const QList<QTreeWidgetItem*> selected = treeWidget->selectedItems();
    foreach (QTreeWidgetItem* item, selected) {
        if (item->type() == ENZYME_ITEM_TYPE) {
            auto enzItem = static_cast<EnzymeItem*>(item);
            selectedAnnotation = enzItem->getEnzymeAnnotation();
            if (selectedAnnotation == nullptr) {
                continue;
            }
            AnnotationSelection* sel = ctx->getAnnotationsSelection();
            sel->clear();
            sel->add(selectedAnnotation);
        }
    }
    if (selectedAnnotation == nullptr) {
        return;
    }

    // Same path as clicking the site on the circular map: scroll every view onto it and
    // select the site's sequence so the circular map, linear map, and sequence view all highlight it.
    const bool crossesJunction = U1AnnotationUtils::isAnnotationContainsJunctionPoint(selectedAnnotation, ctx->getSequenceLength());
    const int regionIndex = crossesJunction ? -1 : 0;
    ctx->emitAnnotationActivated(selectedAnnotation, regionIndex);
    ctx->getSequenceSelection()->setSelectedRegions(selectedAnnotation->getRegions());
}

void RestrctionMapWidget::sl_onAnnotationsGroupCreated(AnnotationGroup* g) {
    if (g->getName() == ANNOTATION_GROUP_ENZYME) {
        updateTreeWidget();
    }
}

void RestrctionMapWidget::initTreeWidget() {
    QSet<AnnotationTableObject*> aObjs = ctx->getAnnotationObjects(true);
    foreach (AnnotationTableObject* obj, aObjs) {
        queueAnnotations(obj->getAnnotations());
    }
}

void RestrctionMapWidget::sl_onAnnotationsInGroupRemoved(const QList<Annotation*>& anns, AnnotationGroup* group) {
    if (group->getName() == ANNOTATION_GROUP_ENZYME) {
        removeAnnotationsFromMap(anns);
    }
}

void RestrctionMapWidget::queueAnnotations(const QList<Annotation*>& anns) {
    if (methylation == nullptr || anns.isEmpty()) {
        return;
    }
    const int queuedBefore = methylation->pendingAnnotations.size();
    for (Annotation* annotation : anns) {
        if (annotation == nullptr) {
            continue;
        }
        if (methylation->itemsByAnnotation.contains(annotation) || methylation->pendingAnnotationLive.contains(annotation) ||
            methylation->heldTokenByAnnotation.contains(annotation)) {
            continue;
        }
        if (!methylation->foldersByName.contains(annotation->getName())) {
            continue;
        }
        methylation->pendingAnnotations.append(annotation);
        methylation->pendingAnnotationLive.insert(annotation);
    }
    if (methylation->pendingAnnotations.size() == queuedBefore || methylation->insertRunning || methylation->insertScheduled) {
        return;
    }
    const int queued = methylation->pendingAnnotations.size() - methylation->pendingAnnotationCursor;
    if (queued > INTERACTIVE_TREE_LIMIT) {
        scheduleAnnotationInsert();
    } else {
        insertAnnotationChunk();
    }
}

void RestrctionMapWidget::scheduleAnnotationInsert() {
    if (methylation == nullptr || methylation->insertScheduled) {
        return;
    }
    methylation->insertScheduled = true;
    const int epoch = methylation->epoch;
    QTimer::singleShot(0, this, [this, epoch]() {
        if (methylation == nullptr || methylation->epoch != epoch) {
            return;
        }
        insertAnnotationChunk();
    });
}

void RestrctionMapWidget::insertAnnotationChunk() {
    if (methylation == nullptr) {
        return;
    }
    methylation->insertScheduled = false;
    if (methylation->pendingAnnotationCursor >= methylation->pendingAnnotations.size()) {
        return;
    }
    methylation->insertRunning = true;
    const bool slice = methylation->pendingAnnotations.size() - methylation->pendingAnnotationCursor > INTERACTIVE_TREE_LIMIT;
    // Turning updates back on lays the whole tree out again. On a long fill that layout
    // outgrows the slice and the cursor stays a pinwheel, so updates stay off until the queue drains.
    if (slice) {
        methylation->bulkFill = true;
        if (!methylation->columnLocked) {
            treeWidget->header()->resizeSection(1, METHYLATION_COLUMN_WIDTH);
            treeWidget->header()->setSectionResizeMode(1, QHeaderView::Fixed);
            methylation->columnLocked = true;
        }
    }
    QElapsedTimer timer;
    timer.start();
    QHash<EnzymeFolderItem*, QList<QTreeWidgetItem*>> batches;
    QSet<EnzymeFolderItem*> heldFolders;
    treeWidget->setUpdatesEnabled(false);
    const bool signalsBlocked = treeWidget->blockSignals(true);
    bool revealedSite = false;
    while (methylation->pendingAnnotationCursor < methylation->pendingAnnotations.size()) {
        if (slice && timer.elapsed() >= INTERACTIVE_SLICE_MS && (!batches.isEmpty() || !heldFolders.isEmpty())) {
            break;
        }
        Annotation* annotation = methylation->pendingAnnotations.at(methylation->pendingAnnotationCursor++);
        if (!methylation->pendingAnnotationLive.contains(annotation)) {
            continue;
        }
        methylation->pendingAnnotationLive.remove(annotation);
        const QString enzymeName = annotation->getName();
        EnzymeFolderItem* folder = methylation->foldersByName.value(enzymeName);
        if (folder == nullptr) {
            continue;
        }
        if (folder->siteCount() == 0 && !batches.contains(folder) && !heldFolders.contains(folder)) {
            revealedSite = true;
        }
        const quint64 token = methylation->nextToken++;
        if (enzymeMayBeBlockedByHostMethylation(enzymeName)) {
            MethylationState::Site site;
            site.token = token;
            site.enzymeName = enzymeName;
            site.regions = annotation->getRegions();
            methylation->pendingSites.append(site);
        }
        int stagedChildren = 0;
        const auto batchIt = batches.constFind(folder);
        if (batchIt != batches.cend()) {
            stagedChildren = batchIt.value().size();
        }
        const int stagedCount = folder->siteCount() + stagedChildren;
        if (stagedCount >= MATERIALIZED_SITES_PER_ENZYME) {
            MethylationState::HeldSite held;
            held.annotation = annotation;
            methylation->heldByToken.insert(token, held);
            methylation->heldTokenByAnnotation.insert(annotation, token);
            methylation->heldTokensByEnzyme[enzymeName].append(token);
            folder->addHeldSite();
            heldFolders.insert(folder);
            continue;
        }
        EnzymeItem* item = folder->createEnzymeItem(annotation, token);
        methylation->itemsByAnnotation.insert(annotation, item);
        methylation->itemsByToken.insert(token, item);
        batches[folder].append(item);
    }
    for (auto it = batches.cbegin(); it != batches.cend(); ++it) {
        EnzymeFolderItem* folder = it.key();
        folder->addChildren(it.value());
        folder->publishSiteCount();
    }
    for (EnzymeFolderItem* folder : heldFolders) {
        if (!batches.contains(folder)) {
            folder->publishSiteCount();
        }
    }
    const bool drained = methylation->pendingAnnotationCursor >= methylation->pendingAnnotations.size();
    if (drained) {
        methylation->pendingAnnotations.clear();
        methylation->pendingAnnotationCursor = 0;
        methylation->pendingAnnotationLive.clear();
        methylation->bulkFill = false;
        // Enzyme rows stay alphabetical. Sites stay in the order they were found.
        // sortChildren does not recurse, so site rows keep discovery order.
        treeWidget->invisibleRootItem()->sortChildren(0, Qt::AscendingOrder);
    }
    treeWidget->blockSignals(signalsBlocked);
    if (!methylation->bulkFill) {
        if (methylation->pendingMarks.isEmpty()) {
            releaseMethylationColumn();
        }
        treeWidget->setUpdatesEnabled(true);
        if (revealedSite || drained) {
            updateEmptyEnzymeVisibility();
        }
    }
    methylation->insertRunning = false;
    if (!drained) {
        scheduleAnnotationInsert();
    }
    if (methylation->pendingSites.isEmpty() || methylation->jobActive) {
        return;
    }
    // A short list is marked before this call returns. A long list is marked on a worker
    // while later slices keep adding rows.
    if (drained && methylation->pendingSites.size() <= INTERACTIVE_TREE_LIMIT) {
        markQueuedSitesOnGuiThread();
    } else if (drained || methylation->pendingSites.size() > INTERACTIVE_TREE_LIMIT) {
        startMethylationJob();
    }
}

void RestrctionMapWidget::removeAnnotationsFromMap(const QList<Annotation*>& anns) {
    if (ctx == nullptr || methylation == nullptr || anns.isEmpty()) {
        return;
    }
    AnnotationSelection* selection = ctx->getAnnotationsSelection();
    QVector<EnzymeItem*> removed;
    removed.reserve(anns.size());
    QSet<EnzymeFolderItem*> heldFolders;
    for (Annotation* annotation : anns) {
        if (annotation == nullptr) {
            continue;
        }
        methylation->pendingAnnotationLive.remove(annotation);
        EnzymeItem* item = methylation->itemsByAnnotation.take(annotation);
        const quint64 heldToken = item == nullptr ? methylation->heldTokenByAnnotation.take(annotation) : 0;
        const bool listedEnzyme = item != nullptr || heldToken != 0 || methylation->foldersByName.contains(annotation->getName());
        if (!listedEnzyme) {
            continue;
        }
        selection->remove(annotation);
        if (item == nullptr) {
            if (heldToken == 0) {
                continue;
            }
            const MethylationState::HeldSite held = methylation->heldByToken.take(heldToken);
            EnzymeFolderItem* folder = methylation->foldersByName.value(annotation->getName());
            if (folder == nullptr) {
                continue;
            }
            folder->adjustBlockedCounts(held.flags, HostMethylationNone);
            folder->removeHeldSite();
            heldFolders.insert(folder);
            continue;
        }
        item->clearAnnotation();
        methylation->itemsByToken.remove(item->getToken());
        removed.append(item);
    }
    if (!heldFolders.isEmpty()) {
        treeWidget->setUpdatesEnabled(false);
        for (EnzymeFolderItem* folder : heldFolders) {
            folder->publishSiteCount();
            folder->publishMethylationMarker();
        }
        if (!methylation->bulkFill) {
            treeWidget->setUpdatesEnabled(true);
        }
    }
    if (removed.isEmpty()) {
        if (!heldFolders.isEmpty() && !methylation->bulkFill) {
            updateEmptyEnzymeVisibility();
        }
        return;
    }
    if (removed.size() <= INTERACTIVE_TREE_LIMIT && !methylation->detachRunning && methylation->pendingDetach.isEmpty()) {
        detachEnzymeItems(removed);
        return;
    }
    methylation->pendingDetach += removed;
    scheduleDetach();
}

void RestrctionMapWidget::detachEnzymeItems(const QVector<EnzymeItem*>& items) {
    if (items.isEmpty() || treeWidget == nullptr) {
        return;
    }
    QSet<EnzymeFolderItem*> dirtyFolders;
    treeWidget->setUpdatesEnabled(false);
    const bool signalsBlocked = treeWidget->blockSignals(true);
    for (EnzymeItem* item : items) {
        QTreeWidgetItem* parent = item->parent();
        if (parent == nullptr || parent->type() != ENZYME_FOLDER_ITEM_TYPE) {
            delete item;
            continue;
        }
        auto* folder = static_cast<EnzymeFolderItem*>(parent);
        folder->detachEnzymeItem(item);
        dirtyFolders.insert(folder);
    }
    for (EnzymeFolderItem* folder : dirtyFolders) {
        folder->publishSiteCount();
        folder->publishMethylationMarker();
    }
    treeWidget->blockSignals(signalsBlocked);
    if (methylation == nullptr || !methylation->bulkFill) {
        treeWidget->setUpdatesEnabled(true);
        updateEmptyEnzymeVisibility();
    }
}

void RestrctionMapWidget::scheduleDetach() {
    if (methylation == nullptr || methylation->detachScheduled) {
        return;
    }
    methylation->detachScheduled = true;
    const int epoch = methylation->epoch;
    QTimer::singleShot(0, this, [this, epoch]() {
        if (methylation == nullptr || methylation->epoch != epoch) {
            return;
        }
        detachNextEnzymeItems();
    });
}

void RestrctionMapWidget::detachNextEnzymeItems() {
    if (methylation == nullptr || treeWidget == nullptr) {
        return;
    }
    methylation->detachScheduled = false;
    if (methylation->pendingDetach.isEmpty()) {
        return;
    }
    methylation->detachRunning = true;
    QSet<EnzymeFolderItem*> dirtyFolders;
    treeWidget->setUpdatesEnabled(false);
    const bool signalsBlocked = treeWidget->blockSignals(true);
    QElapsedTimer timer;
    timer.start();
    while (!methylation->pendingDetach.isEmpty()) {
        EnzymeItem* item = methylation->pendingDetach.takeLast();
        QTreeWidgetItem* parent = item->parent();
        if (parent == nullptr || parent->type() != ENZYME_FOLDER_ITEM_TYPE) {
            delete item;
        } else {
            auto* folder = static_cast<EnzymeFolderItem*>(parent);
            folder->detachEnzymeItem(item);
            dirtyFolders.insert(folder);
        }
        if (timer.elapsed() >= INTERACTIVE_SLICE_MS) {
            break;
        }
    }
    for (EnzymeFolderItem* folder : dirtyFolders) {
        folder->publishSiteCount();
        folder->publishMethylationMarker();
    }
    treeWidget->blockSignals(signalsBlocked);
    if (!methylation->bulkFill) {
        treeWidget->setUpdatesEnabled(true);
        updateEmptyEnzymeVisibility();
    }
    methylation->detachRunning = false;
    if (!methylation->pendingDetach.isEmpty()) {
        scheduleDetach();
    }
}

void RestrctionMapWidget::markQueuedSitesOnGuiThread() {
    if (methylation == nullptr || methylation->pendingSites.isEmpty()) {
        return;
    }
    QVector<MethylationState::Mark> marks;
    if (!methylation->sequenceReady) {
        methylation->sequenceReady = loadSequence(methylation->sequence, methylation->circular);
    }
    if (methylation->sequenceReady && !methylation->sequence.isEmpty()) {
        marks.reserve(methylation->pendingSites.size());
        for (const MethylationState::Site& site : methylation->pendingSites) {
            const int flags = hostMethylationBlock(site.enzymeName, methylation->sequence, methylation->circular, site.regions);
            if (flags == HostMethylationNone) {
                continue;
            }
            MethylationState::Mark mark;
            mark.token = site.token;
            mark.flags = flags;
            marks.append(mark);
        }
    }
    methylation->pendingSites.clear();
    if (marks.isEmpty()) {
        return;
    }
    methylation->pendingMarks += marks;
    applyMethylationChunk();
}

class HostMethylationRunnable : public QRunnable {
public:
    using Mailbox = RestrctionMapWidget::MethylationState::Mailbox;
    using Site = RestrctionMapWidget::MethylationState::Site;
    using Mark = RestrctionMapWidget::MethylationState::Mark;

    HostMethylationRunnable(std::shared_ptr<Mailbox> mailbox, QByteArray sequence, bool circular, QVector<Site> sites)
        : mailbox(std::move(mailbox)), sequence(std::move(sequence)), circular(circular), sites(std::move(sites)) {
        setAutoDelete(true);
    }

    void run() override {
        QVector<Mark> marks;
        marks.reserve(sites.size());
        for (const Site& site : sites) {
            if (mailbox->cancelled.load(std::memory_order_relaxed)) {
                return;
            }
            const int flags = hostMethylationBlock(site.enzymeName, sequence, circular, site.regions);
            if (flags == HostMethylationNone) {
                continue;
            }
            Mark mark;
            mark.token = site.token;
            mark.flags = flags;
            marks.append(mark);
        }
        QMutexLocker locker(&mailbox->mutex);
        if (mailbox->cancelled.load(std::memory_order_relaxed)) {
            return;
        }
        mailbox->marks = std::move(marks);
        mailbox->ready = true;
    }

private:
    std::shared_ptr<Mailbox> mailbox;
    QByteArray sequence;
    bool circular = false;
    QVector<Site> sites;
};

void RestrctionMapWidget::startMethylationJob() {
    if (methylation == nullptr || methylation->jobActive || methylation->pendingSites.isEmpty()) {
        return;
    }
    if (!methylation->sequenceReady) {
        methylation->sequenceReady = loadSequence(methylation->sequence, methylation->circular);
    }
    if (!methylation->sequenceReady || methylation->sequence.isEmpty()) {
        methylation->pendingSites.clear();
        return;
    }
    // The match reads a copy of the sequence off the GUI thread. Rows are already in the
    // tree; the red X labels are applied in short slices after the worker finishes.
    methylation->jobActive = true;
    methylation->mailbox = std::make_shared<MethylationState::Mailbox>();
    auto* runnable = new HostMethylationRunnable(methylation->mailbox, methylation->sequence, methylation->circular, std::move(methylation->pendingSites));
    QThreadPool::globalInstance()->start(runnable);
    scheduleMethylationPoll();
}

void RestrctionMapWidget::scheduleMethylationPoll() {
    if (methylation == nullptr || methylation->pollScheduled) {
        return;
    }
    methylation->pollScheduled = true;
    const int epoch = methylation->epoch;
    QTimer::singleShot(30, this, [this, epoch]() {
        if (methylation == nullptr || methylation->epoch != epoch) {
            return;
        }
        pollMethylation();
    });
}

void RestrctionMapWidget::pollMethylation() {
    if (methylation == nullptr) {
        return;
    }
    methylation->pollScheduled = false;
    if (methylation->mailbox == nullptr) {
        methylation->jobActive = false;
        if (!methylation->pendingSites.isEmpty()) {
            startMethylationJob();
        }
        return;
    }
    QVector<MethylationState::Mark> marks;
    bool ready = false;
    {
        QMutexLocker locker(&methylation->mailbox->mutex);
        if (methylation->mailbox->ready) {
            marks = std::move(methylation->mailbox->marks);
            ready = true;
        }
    }
    if (!ready) {
        scheduleMethylationPoll();
        return;
    }
    methylation->mailbox.reset();
    methylation->jobActive = false;
    if (!marks.isEmpty()) {
        methylation->pendingMarks += marks;
        scheduleMarkApply();
    }
    if (!methylation->pendingSites.isEmpty()) {
        startMethylationJob();
    }
}

void RestrctionMapWidget::scheduleMarkApply() {
    if (methylation == nullptr || methylation->applyScheduled) {
        return;
    }
    methylation->applyScheduled = true;
    const int epoch = methylation->epoch;
    QTimer::singleShot(0, this, [this, epoch]() {
        if (methylation == nullptr || methylation->epoch != epoch) {
            return;
        }
        applyMethylationChunk();
    });
}

void RestrctionMapWidget::applyMethylationChunk() {
    if (methylation == nullptr || treeWidget == nullptr) {
        return;
    }
    methylation->applyScheduled = false;
    if (methylation->pendingMarkCursor >= methylation->pendingMarks.size()) {
        methylation->pendingMarks.clear();
        methylation->pendingMarkCursor = 0;
        releaseMethylationColumn();
        return;
    }
    const bool slice = methylation->pendingMarks.size() - methylation->pendingMarkCursor > INTERACTIVE_TREE_LIMIT;
    if (slice && !methylation->columnLocked) {
        treeWidget->header()->resizeSection(1, METHYLATION_COLUMN_WIDTH);
        treeWidget->header()->setSectionResizeMode(1, QHeaderView::Fixed);
        methylation->columnLocked = true;
    }
    QElapsedTimer timer;
    timer.start();
    QSet<EnzymeFolderItem*> dirtyFolders;
    treeWidget->setUpdatesEnabled(false);
    while (methylation->pendingMarkCursor < methylation->pendingMarks.size()) {
        if (slice && timer.elapsed() >= INTERACTIVE_SLICE_MS) {
            break;
        }
        const MethylationState::Mark& mark = methylation->pendingMarks.at(methylation->pendingMarkCursor++);
        EnzymeItem* item = methylation->itemsByToken.value(mark.token, nullptr);
        if (item == nullptr) {
            auto heldIt = methylation->heldByToken.find(mark.token);
            if (heldIt == methylation->heldByToken.end() || heldIt.value().annotation == nullptr) {
                continue;
            }
            EnzymeFolderItem* folder = methylation->foldersByName.value(heldIt.value().annotation->getName());
            if (folder == nullptr) {
                continue;
            }
            const int oldFlags = heldIt.value().flags;
            if (oldFlags == mark.flags) {
                continue;
            }
            heldIt.value().flags = mark.flags;
            folder->adjustBlockedCounts(oldFlags, mark.flags);
            dirtyFolders.insert(folder);
            continue;
        }
        const int oldFlags = item->getMethylationFlags();
        if (oldFlags == mark.flags) {
            continue;
        }
        QTreeWidgetItem* parent = item->parent();
        if (parent == nullptr || parent->type() != ENZYME_FOLDER_ITEM_TYPE) {
            continue;
        }
        item->setMethylationFlags(mark.flags);
        auto* folder = static_cast<EnzymeFolderItem*>(parent);
        folder->adjustBlockedCounts(oldFlags, mark.flags);
        dirtyFolders.insert(folder);
    }
    for (EnzymeFolderItem* folder : dirtyFolders) {
        folder->publishMethylationMarker();
    }
    restoreTreeUpdates();
    const bool drained = methylation->pendingMarkCursor >= methylation->pendingMarks.size();
    if (drained) {
        methylation->pendingMarks.clear();
        methylation->pendingMarkCursor = 0;
        releaseMethylationColumn();
        return;
    }
    scheduleMarkApply();
}

void RestrctionMapWidget::releaseMethylationColumn() {
    if (methylation == nullptr || !methylation->columnLocked || methylation->bulkFill) {
        return;
    }
    methylation->columnLocked = false;
    if (treeWidget != nullptr) {
        treeWidget->header()->setSectionResizeMode(1, QHeaderView::ResizeToContents);
    }
}

void RestrctionMapWidget::invalidateMethylationTracking() {
    if (methylation == nullptr) {
        return;
    }
    ++methylation->epoch;
    if (methylation->mailbox != nullptr) {
        methylation->mailbox->cancelled.store(true, std::memory_order_relaxed);
        methylation->mailbox.reset();
    }
    methylation->insertRunning = false;
    methylation->insertScheduled = false;
    methylation->detachRunning = false;
    methylation->detachScheduled = false;
    methylation->jobActive = false;
    methylation->pollScheduled = false;
    methylation->applyScheduled = false;
    methylation->pendingAnnotations.clear();
    methylation->pendingAnnotationCursor = 0;
    methylation->pendingAnnotationLive.clear();
    methylation->itemsByAnnotation.clear();
    methylation->itemsByToken.clear();
    methylation->foldersByName.clear();
    methylation->pendingDetach.clear();
    methylation->pendingSites.clear();
    methylation->pendingMarks.clear();
    methylation->pendingMarkCursor = 0;
    methylation->sequenceReady = false;
    methylation->sequence.clear();
    methylation->circular = false;
    methylation->bulkFill = false;
    methylation->materializeScheduled = false;
    methylation->heldByToken.clear();
    methylation->heldTokenByAnnotation.clear();
    methylation->heldTokensByEnzyme.clear();
    methylation->pendingMaterialize.clear();
    methylation->pendingMaterializeCursor = 0;
    releaseMethylationColumn();
    if (treeWidget != nullptr) {
        treeWidget->setUpdatesEnabled(true);
    }
}

void RestrctionMapWidget::restoreTreeUpdates() {
    if (methylation == nullptr || treeWidget == nullptr || methylation->bulkFill) {
        return;
    }
    treeWidget->setUpdatesEnabled(true);
}

void RestrctionMapWidget::sl_itemExpanded(QTreeWidgetItem* item) {
    if (methylation == nullptr || item == nullptr || item->type() != ENZYME_FOLDER_ITEM_TYPE) {
        return;
    }
    auto* folder = static_cast<EnzymeFolderItem*>(item);
    const QVector<quint64> tokens = methylation->heldTokensByEnzyme.take(folder->getName());
    if (tokens.isEmpty()) {
        return;
    }
    methylation->pendingMaterialize += tokens;
    scheduleMaterialize();
}

void RestrctionMapWidget::scheduleMaterialize() {
    if (methylation == nullptr || methylation->materializeScheduled) {
        return;
    }
    methylation->materializeScheduled = true;
    const int epoch = methylation->epoch;
    QTimer::singleShot(0, this, [this, epoch]() {
        if (methylation == nullptr || methylation->epoch != epoch) {
            return;
        }
        materializeHeldChunk();
    });
}

void RestrctionMapWidget::materializeHeldChunk() {
    if (methylation == nullptr || treeWidget == nullptr) {
        return;
    }
    methylation->materializeScheduled = false;
    if (methylation->pendingMaterializeCursor >= methylation->pendingMaterialize.size()) {
        methylation->pendingMaterialize.clear();
        methylation->pendingMaterializeCursor = 0;
        restoreTreeUpdates();
        return;
    }
    const bool slice = methylation->pendingMaterialize.size() - methylation->pendingMaterializeCursor > INTERACTIVE_TREE_LIMIT;
    QElapsedTimer timer;
    timer.start();
    QHash<EnzymeFolderItem*, QList<QTreeWidgetItem*>> batches;
    QSet<EnzymeFolderItem*> touched;
    treeWidget->setUpdatesEnabled(false);
    const bool signalsBlocked = treeWidget->blockSignals(true);
    while (methylation->pendingMaterializeCursor < methylation->pendingMaterialize.size()) {
        if (slice && timer.elapsed() >= INTERACTIVE_SLICE_MS && !batches.isEmpty()) {
            break;
        }
        const quint64 token = methylation->pendingMaterialize.at(methylation->pendingMaterializeCursor++);
        auto heldIt = methylation->heldByToken.find(token);
        if (heldIt == methylation->heldByToken.end()) {
            continue;
        }
        Annotation* annotation = heldIt.value().annotation;
        const int flags = heldIt.value().flags;
        methylation->heldByToken.erase(heldIt);
        if (annotation == nullptr) {
            continue;
        }
        methylation->heldTokenByAnnotation.remove(annotation);
        EnzymeFolderItem* folder = methylation->foldersByName.value(annotation->getName());
        if (folder == nullptr) {
            continue;
        }
        EnzymeItem* item = folder->createEnzymeItem(annotation, token);
        if (flags != 0) {
            // Folder counts were updated when the mark landed on the held site.
            item->setMethylationFlags(flags);
        }
        methylation->itemsByAnnotation.insert(annotation, item);
        methylation->itemsByToken.insert(token, item);
        folder->removeHeldSite();
        batches[folder].append(item);
        touched.insert(folder);
    }
    for (auto it = batches.cbegin(); it != batches.cend(); ++it) {
        it.key()->addChildren(it.value());
    }
    for (EnzymeFolderItem* folder : touched) {
        folder->publishSiteCount();
    }
    treeWidget->blockSignals(signalsBlocked);
    const bool drained = methylation->pendingMaterializeCursor >= methylation->pendingMaterialize.size();
    if (drained) {
        methylation->pendingMaterialize.clear();
        methylation->pendingMaterializeCursor = 0;
    }
    restoreTreeUpdates();
    if (!drained) {
        scheduleMaterialize();
    }
}

}  // namespace U2
