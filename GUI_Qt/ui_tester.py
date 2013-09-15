# -*- coding: utf-8 -*- hwlogic.py
"""
@author: Afsar
This is the Implementation logic for Hello World! (hw.py)
"""
from PySide import QtCore,QtGui
import sys
import qsUI

from enable.api import Window

from chaco.api import ArrayPlotData, Plot
from chaco.tools.api import PanTool, ZoomTool
from numpy import array, linspace
from scipy.special import jn

from pyface.util.guisupport import get_app_qt4, start_event_loop_qt4

class Plotter():
    def __init__(self, parent):
        self.plotdata = ArrayPlotData(x=array([]),  y=array([]))
        self.window = self.create_plot(parent)
        self.widget = self.window.control
        
    def update_data(self, x, y):
        self.plotdata.set_data("x", x)
        self.plotdata.set_data("y", y)

    def create_plot(self, parent):
		x = linspace(-2.0, 10.0, 100)
		pd = ArrayPlotData(index = x)
		for i in range(5):
			pd.set_data("y" + str(i), jn(i,x))

		# Create some line plots of some of the data
		plot = Plot(pd, padding=[40,10,0,40], border_visible=True)
		plot.legend.visible = True
		plot.plot(("index", "y0", "y1", "y2"), name="j_n, n<3", color="red")
		plot.plot(("index", "y3"), name="j_3", color="blue")

		# Attach some tools to the plot
		plot.tools.append(PanTool(plot))
		zoom = ZoomTool(component=plot, tool_mode="box", always_on=False)
		plot.overlays.append(zoom)

		# This Window object bridges the Enable and Qt4 worlds, and handles events
		# and drawing.  We can create whatever hierarchy of nested containers we
		# want, as long as the top-level item gets set as the .component attribute
		# of a Window.
		return Window(parent, -1, component=plot)

class PlotWindow(QtGui.QMainWindow,qsUI.Ui_MainWindow):
    """
    hwl is inherited from both QtGui.QDialog and hw.Ui_Dialog
    """
    def __init__(self,parent=None):
		"""
		Initialization of the class. Call the __init__ for the super classes
		"""
		super(PlotWindow,self).__init__(parent)
		self.setupUi(self)
		self.plotter=Plotter(self)
		self.plotLayout.addWidget(self.plotter.widget)
		self.plotWidget=self.plotter.widget
		self.plotWidget.show()
        #self.connectActions()
		
    def main(self):
        self.show()
	
	'''
    def connectActions(self):
        """
        Connect the user interface controls to the logic
        """
        self.cmdWrite.clicked.connect(self.myprint)      
        
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
	hwl1 = PlotWindow()
	hwl1.main()
	start_event_loop_qt4(app)
	#sys.exit(app.exec_())