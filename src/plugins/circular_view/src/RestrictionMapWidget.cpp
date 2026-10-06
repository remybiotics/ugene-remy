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

#include <QButtonGroup>
#include <QHBoxLayout>
#include <QPushButton>
#include <QSizePolicy>
#include <QVBoxLayout>

#include <U2Algorithm/EnzymeModel.h>

#include <U2Core/Annotation.h>
#include <U2Core/AnnotationGroup.h>
#include <U2Core/AnnotationSelection.h>
#include <U2Core/AnnotationTableObject.h>
#include <U2Core/AppContext.h>
#include <U2Core/AutoAnnotationsSupport.h>
#include <U2Core/DNAAlphabet.h>
#include <U2Core/DNASequenceSelection.h>
#include <U2Core/Settings.h>
#include <U2Core/U1AnnotationUtils.h>

#include <U2Gui/GUIUtils.h>

#include <U2View/ADVSequenceObjectContext.h>
#include <U2View/AnnotatedDNAView.h>
#include <U2View/AutoAnnotationUtils.h>

#define ENZYME_FOLDER_ITEM_TYPE 1022
#define ENZYME_ITEM_TYPE 1023

namespace {

const char* SHOW_EMPTY_ENZYMES_KEY = "circular_view/restriction_map_show_empty_enzymes";
const int UNLIMITED_HIT_COUNT = INT_MAX;

}  // namespace

namespace U2 {

//////////////////////////////////////////////////////////////////////////
/// EnzymeItem

EnzymeItem::EnzymeItem(const QString& location, Annotation* a)
    : QTreeWidgetItem(QStringList(location), ENZYME_ITEM_TYPE), annotation(a) {
}

//////////////////////////////////////////////////////////////////////////
/// EnzymeFolderItem

EnzymeFolderItem::EnzymeFolderItem(const QString& name)
    : QTreeWidgetItem(ENZYME_FOLDER_ITEM_TYPE), enzymeName(name) {
    setText(0, QString("%1 : %2 %3").arg(name).arg(0).arg("sites"));
}

void EnzymeFolderItem::addEnzymeItem(Annotation* enzAnn) {
    const SharedAnnotationData& data = enzAnn->getData();
    QString location = U1AnnotationUtils::buildLocationString(data);
    addChild(new EnzymeItem(location, enzAnn));
    setIcon(0, QIcon(":circular_view/images/folder.png"));
    int count = childCount();
    QString site = count == 1 ? RestrctionMapWidget::tr("site") : RestrctionMapWidget::tr("sites");
    setText(0, QString("%1 : %2 %3").arg(getName()).arg(count).arg(site));
}

void EnzymeFolderItem::removeEnzymeItem(Annotation* enzAnn) {
    int count = childCount();
    for (int i = 0; i < count; ++i) {
        auto item = static_cast<EnzymeItem*>(child(i));
        if (item->getEnzymeAnnotation() == enzAnn) {
            removeChild(item);
            QString site = --count == 1 ? RestrctionMapWidget::tr("site") : RestrctionMapWidget::tr("sites");
            setText(0, QString("%1 : %2 %3").arg(getName()).arg(count).arg(site));
            if (count == 0) {
                setIcon(0, QIcon(":circular_view/images/empty_folder.png"));
            }
            break;
        }
    }
}

//////////////////////////////////////////////////////////////////////////
/// RestrictionMapWidget

RestrctionMapWidget::RestrctionMapWidget(ADVSequenceObjectContext* context, QWidget* p)
    : QWidget(p), ctx(nullptr), annotatedDnaView(nullptr), treeWidget(nullptr) {
    assert(context != nullptr);
    setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Expanding);
    annotatedDnaView = context->getAnnotatedDNAView();

    auto layout = new QVBoxLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);

    treeWidget = new QTreeWidget(this);
    treeWidget->setObjectName("restrictionMapTreeWidget");
    treeWidget->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Expanding);
    treeWidget->setColumnCount(1);
    treeWidget->setHeaderLabel(tr("Restriction Sites Map"));
    connect(treeWidget, SIGNAL(itemSelectionChanged()), SLOT(sl_itemSelectionChanged()));

    layout->addWidget(treeWidget, 1);
    initFooterButtons();

    connect(annotatedDnaView, SIGNAL(si_activeSequenceWidgetChanged(ADVSequenceWidget*, ADVSequenceWidget*)), SLOT(sl_onActiveSequenceChanged()));
    connect(annotatedDnaView, SIGNAL(si_sequenceRemoved(ADVSequenceObjectContext*)), SLOT(sl_onSequenceRemoved(ADVSequenceObjectContext*)));
    setSequenceContext(context);
}

void RestrctionMapWidget::setSequenceContext(ADVSequenceObjectContext* sequenceContext) {
    if (ctx == sequenceContext) {
        return;
    }
    unregisterAnnotationObjects();
    ctx = sequenceContext;
    if (ctx == nullptr) {
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
        const bool empty = item->childCount() == 0;
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

void RestrctionMapWidget::sl_onAnnotationsAdded(const QList<Annotation*>& anns) {
    foreach (Annotation* a, anns) {
        QString aName = a->getName();
        EnzymeFolderItem* folderItem = findEnzymeFolderByName(aName);
        if (folderItem) {
            folderItem->addEnzymeItem(a);
        }
    }

    // TODO: enable "intelligent" sorting by reimplementing custom AbstractModel
    //  Take into account number of items in each enzymes folder
    treeWidget->sortItems(0, Qt::AscendingOrder);
    updateEmptyEnzymeVisibility();
}

void RestrctionMapWidget::sl_onAnnotationsRemoved(const QList<Annotation*>& anns) {
    if (ctx == nullptr) {
        return;
    }
    foreach (Annotation* a, anns) {
        EnzymeFolderItem* folderItem = findEnzymeFolderByName(a->getName());
        if (folderItem) {
            AnnotationSelection* sel = ctx->getAnnotationsSelection();
            sel->remove(a);
            folderItem->removeEnzymeItem(a);
        }
    }
    updateEmptyEnzymeVisibility();
}

EnzymeFolderItem* RestrctionMapWidget::findEnzymeFolderByName(const QString& enzymeName) {
    int count = treeWidget->topLevelItemCount();

    for (int i = 0; i < count; i++) {
        assert(treeWidget->topLevelItem(i)->type() == ENZYME_FOLDER_ITEM_TYPE);
        auto item = static_cast<EnzymeFolderItem*>(treeWidget->topLevelItem(i));
        if (item->getName() == enzymeName) {
            return item;
        }
    }

    return nullptr;
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
        QList<Annotation*> anns = obj->getAnnotations();
        for (Annotation* a : qAsConst(anns)) {
            QString aName = a->getName();
            EnzymeFolderItem* folderItem = findEnzymeFolderByName(aName);
            if (folderItem) {
                folderItem->addEnzymeItem(a);
            }
        }
    }
    treeWidget->sortItems(0, Qt::AscendingOrder);
}

void RestrctionMapWidget::sl_onAnnotationsInGroupRemoved(const QList<Annotation*>& anns, AnnotationGroup* group) {
    if (ctx == nullptr) {
        return;
    }
    if (group->getName() == ANNOTATION_GROUP_ENZYME) {
        foreach (Annotation* a, anns) {
            EnzymeFolderItem* folderItem = findEnzymeFolderByName(a->getName());
            if (folderItem) {
                AnnotationSelection* sel = ctx->getAnnotationsSelection();
                sel->remove(a);
                folderItem->removeEnzymeItem(a);
            }
        }
        updateEmptyEnzymeVisibility();
    }
}

}  // namespace U2
