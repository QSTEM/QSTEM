# -*- coding: utf-8 -*- 
"""
/*
QSTEM - image simulation for TEM/STEM/CBED
    Copyright (C) 2000-2010  Christoph Koch
    Copyright (C) 2010-2013  Christoph Koch, Michael Sarahan

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
"""

from pyface.util.guisupport import get_app_qt4, start_event_loop_qt4

from dialogs.main_window import QSTEM_UI

if __name__ == '__main__':
    #app = QtGui.QApplication(sys.argv)
    app = get_app_qt4()
    gui = QSTEM_UI()
    gui.main()
    start_event_loop_qt4(app)
