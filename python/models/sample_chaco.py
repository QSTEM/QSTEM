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

from chaco.api import ArrayPlotData, Plot
from chaco.tools.api import PanTool, ZoomTool
from traits.api import Instance, on_trait_change
from enable.api import Window

from raw_models import sample

class SamplePlotter(sample.SampleModel):
    pd = Instance(ArrayPlotData)
    plot = Instance(Plot)
    window = Instance(Window)
    
    def __init__(self, parent):
        self.pd = ArrayPlotData()
        self.window = self.create_plot(parent)
        self.widget = self.window.control    
    
    def create_plot(self, parent):
	# Create some line plots of some of the data
	self.plot = Plot(self.pd, padding=[40,10,0,40], border_visible=True)
	self.plot.legend.visible = True
    
	# Attach some tools to the plot
	self.plot.tools.append(PanTool(self.plot))
	zoom = ZoomTool(component=self.plot, tool_mode="box", always_on=False)
	self.plot.overlays.append(zoom)
	
	# This Window object bridges the Enable and Qt4 worlds, and handles events
	# and drawing.  We can create whatever hierarchy of nested containers we
	# want, as long as the top-level item gets set as the .component attribute
	# of a Window.
	return Window(parent, -1, component=self.plot)
    
    def top_plot(self):
	# clear any existing plots
	
	
	#self.plot.aspect_ratio=self.nCellsX
	for key in self.transformed_elements.keys():
	    if key in self.plot.plots:
		self.plot.delplot(key)
	    self.plot.plot(("x"+key,"y"+key), name=key, type="scatter")
	    
    def front_plot(self):
	# clear any existing plots
	for key in self.transformed_elements.keys():
	    if key in self.plot.plots:
		self.plot.delplot(key)	    
	    self.plot.plot(("x"+key,"z"+key), name=key, type="scatter")	
	    
    @on_trait_change("")
    def _update_coordinates(self):
	for key in self.transformed_elements.keys():
	    self.pd.set_data("x" + key, self.transformed_elements[key][:,0])
	    self.pd.set_data("y" + key, self.transformed_elements[key][:,1])
	    self.pd.set_data("z" + key, self.transformed_elements[key][:,2])