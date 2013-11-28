These files are generated using either pyuic4 (provided by PyQt) or pyside-uic (provided by PySide).

They require one manual change whenever you regenerate these files using those tools: lines that are specific to PySide or PyQt should use pyface.qt instead:

Instead of:
from PySide import QtCore,QtGui

Do this:
from pyface.qt import QtCore,QtGui

This keeps the UI generic, and also is important for making the widgets communicate nicely with Chaco, which is used for the plotting.