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

#include "depthmapX/apptheme.h"

#include <QApplication>
#include <QColor>
#include <QFont>
#include <QPalette>
#include <QStyle>
#include <QStyleFactory>

namespace {
    QString s_nativeStyleName;

    struct Colours {
        QColor window, base, alternateBase, text, mutedText, button, border, highlight, highlightText, mdi;
    };

    Colours lightColours() {
        return {QColor("#f5f6f8"), QColor("#ffffff"), QColor("#f0f2f5"), QColor("#1f2329"), QColor("#8a919c"),
                QColor("#ffffff"), QColor("#d5d9e0"), QColor("#2f6fed"), QColor("#ffffff"), QColor("#e4e7ec")};
    }

    Colours darkColours() {
        return {QColor("#1e2127"), QColor("#16181d"), QColor("#23262d"), QColor("#e3e6eb"), QColor("#7d8590"),
                QColor("#2a2e36"), QColor("#3a3f4a"), QColor("#4c8dff"), QColor("#ffffff"), QColor("#121418")};
    }

    QPalette makePalette(const Colours &c) {
        QPalette p;
        p.setColor(QPalette::Window, c.window);
        p.setColor(QPalette::WindowText, c.text);
        p.setColor(QPalette::Base, c.base);
        p.setColor(QPalette::AlternateBase, c.alternateBase);
        p.setColor(QPalette::ToolTipBase, c.base);
        p.setColor(QPalette::ToolTipText, c.text);
        p.setColor(QPalette::Text, c.text);
        p.setColor(QPalette::Button, c.button);
        p.setColor(QPalette::ButtonText, c.text);
        p.setColor(QPalette::BrightText, QColor("#ff5555"));
        p.setColor(QPalette::Highlight, c.highlight);
        p.setColor(QPalette::HighlightedText, c.highlightText);
        p.setColor(QPalette::Link, c.highlight);
        p.setColor(QPalette::Light, c.button.lighter(115));
        p.setColor(QPalette::Midlight, c.border.lighter(110));
        p.setColor(QPalette::Mid, c.border);
        p.setColor(QPalette::Dark, c.border.darker(130));
        p.setColor(QPalette::Shadow, QColor(0, 0, 0, 90));
        p.setColor(QPalette::Disabled, QPalette::Text, c.mutedText);
        p.setColor(QPalette::Disabled, QPalette::WindowText, c.mutedText);
        p.setColor(QPalette::Disabled, QPalette::ButtonText, c.mutedText);
        p.setColor(QPalette::Disabled, QPalette::Highlight, c.border);
        return p;
    }

    QString makeStyleSheet(const Colours &c) {
        // Kept deliberately small: Fusion + the palette do most of the work, this only adds
        // spacing, flat toolbars and clearer separation between panels.
        QString css = R"(
QMainWindow::separator { background: %border; width: 1px; height: 1px; }
QMenuBar { background: %window; border-bottom: 1px solid %border; padding: 2px 4px; }
QMenuBar::item { padding: 4px 10px; border-radius: 4px; background: transparent; }
QMenuBar::item:selected { background: %alt; }
QMenu { background: %base; border: 1px solid %border; padding: 4px; }
QMenu::item { padding: 5px 28px 5px 24px; border-radius: 4px; }
QMenu::item:selected { background: %highlight; color: %highlightText; }
QMenu::item:disabled { color: %muted; }
QMenu::separator { height: 1px; background: %border; margin: 4px 8px; }
QToolBar { background: %window; border: none; border-bottom: 1px solid %border; padding: 3px; spacing: 2px; }
QToolBar::separator { background: %border; width: 1px; margin: 4px 6px; }
QToolButton { border: 1px solid transparent; border-radius: 4px; padding: 3px; }
QToolButton:hover { background: %alt; border-color: %border; }
QToolButton:checked, QToolButton:pressed { background: %highlightSoft; border-color: %highlight; }
QToolButton[popupMode="2"] { padding-right: 14px; }
QDockWidget::title { background: %window; padding: 6px 8px; border-bottom: 1px solid %border; font-weight: 600; }
QTreeView, QTreeWidget, QListView, QTableView, QTableWidget { background: %base; border: 1px solid %border;
    alternate-background-color: %alt; selection-background-color: %highlight; selection-color: %highlightText; }
QHeaderView::section { background: %window; border: none; border-right: 1px solid %border;
    border-bottom: 1px solid %border; padding: 4px 6px; font-weight: 600; }
QTabBar::tab { background: %window; border: 1px solid %border; border-bottom: none; padding: 5px 12px;
    border-top-left-radius: 4px; border-top-right-radius: 4px; margin-right: 2px; }
QTabBar::tab:selected { background: %base; }
QStatusBar { background: %window; border-top: 1px solid %border; }
QStatusBar::item { border: none; }
QPushButton { background: %button; border: 1px solid %border; border-radius: 4px; padding: 5px 14px; }
QPushButton:hover { border-color: %highlight; }
QPushButton:default { border-color: %highlight; }
QPushButton:disabled { color: %muted; }
QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox { background: %base; border: 1px solid %border;
    border-radius: 4px; padding: 3px 6px; }
QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus { border-color: %highlight; }
QGroupBox { border: 1px solid %border; border-radius: 6px; margin-top: 12px; padding-top: 6px; }
QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 4px; font-weight: 600; }
QProgressBar { border: 1px solid %border; border-radius: 4px; background: %base; text-align: center; }
QProgressBar::chunk { background: %highlight; border-radius: 3px; }
QToolTip { background: %base; color: %text; border: 1px solid %border; padding: 4px; }
)";
        QColor soft = c.highlight;
        soft.setAlpha(60);
        css.replace("%highlightText", c.highlightText.name());
        css.replace("%highlightSoft", QString("rgba(%1,%2,%3,%4)")
                                          .arg(soft.red())
                                          .arg(soft.green())
                                          .arg(soft.blue())
                                          .arg(soft.alpha()));
        css.replace("%highlight", c.highlight.name());
        css.replace("%window", c.window.name());
        css.replace("%border", c.border.name());
        css.replace("%base", c.base.name());
        css.replace("%alt", c.alternateBase.name());
        css.replace("%muted", c.mutedText.name());
        css.replace("%button", c.button.name());
        css.replace("%text", c.text.name());
        return css;
    }
} // namespace

namespace AppTheme {

    QStringList availableThemes() { return QStringList() << light << dark << classic; }

    QString displayName(const QString &theme) {
        if (theme == dark)
            return QApplication::translate("AppTheme", "Dark");
        if (theme == classic)
            return QApplication::translate("AppTheme", "Classic (system)");
        return QApplication::translate("AppTheme", "Light");
    }

    void apply(const QString &theme) {
        auto *app = qobject_cast<QApplication *>(QApplication::instance());
        if (!app) {
            return;
        }
        if (s_nativeStyleName.isEmpty() && app->style()) {
            s_nativeStyleName = app->style()->objectName();
        }

        if (theme == classic) {
            if (!s_nativeStyleName.isEmpty()) {
                app->setStyle(QStyleFactory::create(s_nativeStyleName));
            }
            app->setPalette(app->style()->standardPalette());
            app->setStyleSheet(QString());
            return;
        }

        const Colours colours = (theme == dark) ? darkColours() : lightColours();
        app->setStyle(QStyleFactory::create("Fusion"));
        app->setPalette(makePalette(colours));
        app->setStyleSheet(makeStyleSheet(colours));
    }
} // namespace AppTheme
