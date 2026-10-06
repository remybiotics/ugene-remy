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

#pragma once

#include <QTreeWidget>

class QPushButton;

namespace U2 {

class Annotation;
class AnnotationGroup;
class AnnotationTableObject;
class ADVSequenceObjectContext;
class AnnotatedDNAView;

class EnzymeItem : public QTreeWidgetItem {
public:
    EnzymeItem(const QString& locationStr, Annotation* a);
    Annotation* getEnzymeAnnotation() const {
        return annotation;
    }

private:
    Annotation* annotation;
};

class EnzymeFolderItem : public QTreeWidgetItem {
    QString enzymeName;

public:
    EnzymeFolderItem(const QString& name);
    void addEnzymeItem(Annotation* enzAnn);
    void removeEnzymeItem(Annotation* enzAnn);
    const QString& getName() const {
        return enzymeName;
    }
};

class RestrctionMapWidget : public QWidget {
    Q_OBJECT
public:
    RestrctionMapWidget(ADVSequenceObjectContext* ctx, QWidget* p);

private slots:
    void sl_onAnnotationsAdded(const QList<Annotation*>& anns);
    void sl_onAnnotationsRemoved(const QList<Annotation*>& anns);
    void sl_onAnnotationsInGroupRemoved(const QList<Annotation*>& anns, AnnotationGroup* group);
    void sl_onAnnotationsGroupCreated(AnnotationGroup* g);
    void sl_itemSelectionChanged();
    void sl_onActiveSequenceChanged();
    void sl_onSequenceRemoved(ADVSequenceObjectContext* sequenceContext);
    void sl_onAnnotationObjectAdded(AnnotationTableObject* object);
    void sl_onAnnotationObjectRemoved(AnnotationTableObject* object);
    void sl_maxHitsClicked();
    void sl_showEmptyEnzymesToggled(bool showEmpty);

private:
    void setSequenceContext(ADVSequenceObjectContext* sequenceContext);
    void unregisterAnnotationObjects();
    void connectAnnotationObject(AnnotationTableObject* object);
    void rebuildTree();
    void updateEmptyEnzymeVisibility();
    void initFooterButtons();
    int selectedMaxHitCount() const;
    ADVSequenceObjectContext* ctx;
    AnnotatedDNAView* annotatedDnaView;
    QTreeWidget* treeWidget;
    QPushButton* maxHits1Button = nullptr;
    QPushButton* maxHits2Button = nullptr;
    QPushButton* maxHitsAllButton = nullptr;
    bool showEmptyEnzymes = true;
    EnzymeFolderItem* findEnzymeFolderByName(const QString& enzymeName);
    void registerAnnotationObjects();
    void updateTreeWidget();
    void initTreeWidget();
};

}  // namespace U2
