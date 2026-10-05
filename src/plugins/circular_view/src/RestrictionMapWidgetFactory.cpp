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

#include "RestrictionMapWidgetFactory.h"

#include <U2Core/DNAAlphabet.h>
#include <U2Core/U2SafePoints.h>

#include <U2View/ADVSequenceObjectContext.h>
#include <U2View/AnnotatedDNAView.h>

#include "RestrictionMapWidget.h"

namespace U2 {

const QString RestrictionMapWidgetFactory::GROUP_ID = "OP_RESTRICTION_SITES";
const QString RestrictionMapWidgetFactory::GROUP_ICON_STR = ":circular_view/images/side_list.png";
const QString RestrictionMapWidgetFactory::GROUP_DOC_PAGE = "65929523";

RestrictionMapWidgetFactory::RestrictionMapWidgetFactory() {
    objectViewOfWidget = ObjViewType_SequenceView;
}

QWidget* RestrictionMapWidgetFactory::createWidget(GObjectViewController* objView, const QVariantMap& /*options*/) {
    SAFE_POINT(objView != nullptr,
               QString("Internal error: unable to create widget for group '%1', object view is NULL.").arg(GROUP_ID),
               nullptr);

    auto annotatedDnaView = qobject_cast<AnnotatedDNAView*>(objView);
    SAFE_POINT(annotatedDnaView != nullptr,
               QString("Internal error: unable to cast object view to AnnotatedDNAView for group '%1'.").arg(GROUP_ID),
               nullptr);

    ADVSequenceObjectContext* sequenceContext = annotatedDnaView->getActiveSequenceContext();
    if (sequenceContext == nullptr || sequenceContext->getAlphabet() == nullptr || !sequenceContext->getAlphabet()->isNucleic()) {
        sequenceContext = nullptr;
        const QList<ADVSequenceObjectContext*> contexts = annotatedDnaView->getSequenceContexts();
        foreach (ADVSequenceObjectContext* candidate, contexts) {
            if (candidate->getAlphabet() != nullptr && candidate->getAlphabet()->isNucleic()) {
                sequenceContext = candidate;
                break;
            }
        }
    }
    SAFE_POINT(sequenceContext != nullptr, "Restriction sites map requires a nucleic sequence", nullptr);

    auto widget = new RestrctionMapWidget(sequenceContext, nullptr);
    widget->setObjectName("RestrictionMapWidget");
    return widget;
}

OPGroupParameters RestrictionMapWidgetFactory::getOPGroupParameters() {
    return OPGroupParameters(GROUP_ID, GROUP_ICON_STR, QObject::tr("Restriction Sites Map"), GROUP_DOC_PAGE);
}

bool RestrictionMapWidgetFactory::passFiltration(OPFactoryFilterVisitorInterface* filter) {
    SAFE_POINT(filter != nullptr, "OPFactoryFilterVisitorInterface::filter is NULL", false);
    return filter->typePass(getObjectViewType()) && filter->atLeastOneAlphabetPass(DNAAlphabet_NUCL);
}

const QString& RestrictionMapWidgetFactory::getGroupId() {
    return GROUP_ID;
}

}  // namespace U2
