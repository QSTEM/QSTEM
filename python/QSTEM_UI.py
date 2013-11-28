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
from pyface.qt import QtCore,QtGui
import sys
from dialogs import main_window

from pyface.util.guisupport import get_app_qt4, start_event_loop_qt4

from models.sample_chaco import SamplePlotter

class QSTEM_UI(QtGui.QMainWindow,main_window.Ui_MainWindow):
    """
    hwl is inherited from both QtGui.QDialog and hw.Ui_Dialog
    """
    def __init__(self,parent=None):
        """
        Initialization of the class. Call the __init__ for the super classes
        """
        super(QSTEM_UI,self).__init__(parent)
        self.setupUi(self)
        self.sample=SamplePlotter(self)
        self.plotLayout.addWidget(self.sample.widget)
        self.sampleWidget=self.sample.widget
        self.sampleWidget.show()
        self.connectActions()

    def main(self):
        self.show()

    def connectActions(self):
        """
        Connect the user interface controls to the logic
        """
        self.LoadModel.clicked.connect(self.sample.load_file)

        '''
        def myprint(self):
        """
        Even handler for the pushButton click
        """
        self.txtLine.setText('Python -- ')        
        self.txtEdit.setText('This')
        self.lblShow.setText('is a test')
        '''
        
    
if __name__=='__main__':
    #app = QtGui.QApplication(sys.argv)
    app = get_app_qt4()
    gui = QSTEM_UI()
    gui.main()
    start_event_loop_qt4(app)
    #sys.exit(app.exec_())