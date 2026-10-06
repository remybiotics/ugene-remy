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

#include "GroupOptionsWidget.h"

#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>

#include <U2Gui/HelpButton.h>

namespace U2 {

GroupOptionsWidget::GroupOptionsWidget(const QString& _groupId, const QString& _title, const QString& documentationPage, QWidget* _widget, QWidget* mainWidget)
    : groupId(_groupId),
      widget(_widget),
      mainWidget(mainWidget),
      title(_title) {
#ifdef Q_OS_DARWIN
    setStyleSheet("font-size: 11.25pt;");
#else
    setStyleSheet("font-size: 8.25pt;");
#endif

    titleWidget = new QLabel(title);
    titleWidget->setObjectName("titleWidget");
    titleWidget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
    titleWidget->setMinimumWidth(MIN_WIDGET_WIDTH);

    titleWidget->setStyleSheet(
        "background: palette(midlight);"
        "border-style: solid;"
        "border-width: 1px;"
        "border-color: palette(mid);"
        "padding: 2px;"
        "margin: 5px;");

    widget->setContentsMargins(10, 5, 5, 5);

    // The content takes the leftover height so list tabs fill the options panel the same way Search in Sequence does.
    mainLayout = new QVBoxLayout();
    mainLayout->setContentsMargins(0, 0, 0, 15);
    mainLayout->setSpacing(0);
    mainLayout->addWidget(titleWidget);
    mainLayout->addWidget(widget, 1);

    auto helpButton = new QPushButton(tr("Help"), this);
    helpButton->setMaximumWidth(60);
    new HelpButton(this, helpButton, documentationPage);

    // A content widget can leave a direct child named optionsPanelFooter to sit in the empty space beside Help.
    QWidget* footer = mainWidget == nullptr ? nullptr : mainWidget->findChild<QWidget*>("optionsPanelFooter", Qt::FindDirectChildrenOnly);

    auto helpLayout = new QHBoxLayout();
    helpLayout->setContentsMargins(footer == nullptr ? 0 : 10, 0, 10, 0);
    helpLayout->setSpacing(6);
    if (footer != nullptr) {
        helpLayout->addWidget(footer, 0, Qt::AlignVCenter);
    }
    helpLayout->addStretch(1);
    helpLayout->addWidget(helpButton, 0, Qt::AlignRight | Qt::AlignVCenter);

    mainLayout->addLayout(helpLayout);
    mainLayout->setAlignment(helpLayout, Qt::AlignBottom);

    setLayout(mainLayout);

    setFocusProxy(widget);
}

}  // namespace U2
