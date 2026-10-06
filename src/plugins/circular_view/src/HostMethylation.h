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

#include <QString>
#include <QVector>

#include <U2Core/U2Region.h>

namespace U2 {

enum HostMethylationBlock {
    HostMethylationNone = 0,
    HostMethylationDam = 1,
    HostMethylationDcm = 2
};

/**
 * Dam (GATC) or Dcm (CCWGG) methylation that blocks this enzyme hit in a Dam+ or Dcm+ host.
 * Returns a HostMethylationBlock bitmask. Enzymes that cut methylated DNA, and enzymes that
 * are only slowed by methylation, are not marked.
 *
 * Reads only the arguments and the process-wide rule table, so a worker thread may call it
 * on its own copies of the sequence and regions.
 */
int hostMethylationBlock(const QString& enzymeName, const QByteArray& sequence, bool circular, const QVector<U2Region>& regions);

/** True when this enzyme has a Dam or Dcm block rule. Does not read the sequence. */
bool enzymeMayBeBlockedByHostMethylation(const QString& enzymeName);

}  // namespace U2
