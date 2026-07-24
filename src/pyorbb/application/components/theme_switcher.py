from PySide6 import QtCore

# on macos we can use Cocoa to force dark/light mode
try:
	from AppKit import NSApplication, NSAppearance

	class ThemeSwitcher:
	    LIGHT = "NSAppearanceNameAqua"
	    DARK = "NSAppearanceNameDarkAqua"

	    def __init__(self, parent):
	    	self.parent = parent

	    def set_dark_theme(self):
	        NSApp = NSApplication.sharedApplication()
	        appearance = NSAppearance.appearanceNamed_(self.DARK)
	        NSApp.setAppearance_(appearance)

	    def set_light_theme(self):
	        NSApp = NSApplication.sharedApplication()
	        appearance = NSAppearance.appearanceNamed_(self.LIGHT)
	        NSApp.setAppearance_(appearance)

	    def set_auto_theme(self):
	        NSApp = NSApplication.sharedApplication()
	        appearance = NSAppearance.appearanceNamed_(None)
	        NSApp.setAppearance_(appearance)

except:
	# on other platforms this is not possible yet
	class ThemeSwitcher:
	    def __init__(self, parent):
	    	self.parent = parent

	    def set_dark_theme(self):
	    	...
	        # self.parent.changeEvent(QtCore.QEvent(QtCore.QEvent.Type.ThemeChange))
	        # self.parent.set_style()

	    def set_light_theme(self):
	    	...
	        # self.parent.changeEvent(QtCore.QEvent(QtCore.QEvent.Type.ThemeChange))
	        # self.parent.set_style()

	    def set_auto_theme(self):
	    	...
	        # self.parent.changeEvent(QtCore.QEvent(QtCore.QEvent.Type.ThemeChange))
	        # self.parent.set_style()
