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

#define ENZYME_FOLDER_ITEM_TYPE 1022
#define ENZYME_ITEM_TYPE 1023

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
    }
}

}  // namespace U2
