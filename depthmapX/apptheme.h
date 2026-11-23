// depthmapX - spatial network analysis platform
// Copyright (C) 2026, depthmapX contributors

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <QString>
#include <QStringList>

// Application-wide look and feel. The map views draw their own background (see the
// foreground/background colour options), so the theme only changes the window chrome:
// menus, toolbars, docks, dialogs and tables.
namespace AppTheme {
    const QString light = "light";
    const QString dark = "dark";
    const QString classic = "classic"; // the platform's native style, as in earlier versions

    const QString settingTag = "uiTheme";

    QStringList availableThemes();
    QString displayName(const QString &theme);

    // Applies the theme to the running QApplication. Unknown names fall back to light.
    void apply(const QString &theme);
} // namespace AppTheme
