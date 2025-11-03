import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg
from PySide6 import QtGui, QtCore
from PIL import ImageQt, Image


def convert_to_QPixMap(text, fs=9, darkmode=False, scale=2):

    #---- set up a mpl figure instance ----

    fig = plt.Figure()
    fig.patch.set_facecolor('none')
    fig.set_canvas(FigureCanvasAgg(fig))
    renderer = fig.canvas.get_renderer()

    #---- plot the text expression ----

    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    ax.patch.set_facecolor('none')
    if darkmode:
      t = ax.text(0, 0, text, ha='left', va='bottom', fontsize=fs*scale, c='white')
    else:
      t = ax.text(0, 0, text, ha='left', va='bottom', fontsize=fs*scale)

    #---- fit figure size to text artist ----

    fwidth, fheight = fig.get_size_inches()
    fig_bbox = fig.get_window_extent(renderer)

    text_bbox = t.get_window_extent(renderer)

    tight_fwidth = text_bbox.width * fwidth / fig_bbox.width
    tight_fheight = text_bbox.height * fheight / fig_bbox.height

    fig.set_size_inches(tight_fwidth, tight_fheight)


    #---- convert mpl figure to QPixmap ----

    buf, size = fig.canvas.print_to_buffer()
    qimage = QtGui.QImage.rgbSwapped(QtGui.QImage(buf, size[0], size[1],
                                                  QtGui.QImage.Format_ARGB32))
    qpixmap = QtGui.QPixmap(qimage)
    qpixmap.setDevicePixelRatio(scale)
    return qpixmap


def convert_to_QIcon(text, prefix_icon=None, fs=9, darkmode=False, scale=2):
    qpixmap = convert_to_QPixMap(text, fs=fs, darkmode=darkmode, scale=scale)
    qicon = QtGui.QIcon(qpixmap)

    if prefix_icon is None:
        return qicon

    prefix_icon_image = ImageQt.fromqimage(prefix_icon.pixmap(prefix_icon.actualSize(QtCore.QSize(1024, 1024))).toImage())
    image = ImageQt.fromqimage(qicon.pixmap(qicon.actualSize(QtCore.QSize(1024, 1024))).toImage())


    (width1, height1) = prefix_icon_image.size
    (width2, height2) = image.size

    min_height = min(height1, height2)

    r = height2 / min_height
    image = image.resize((int(image.size[0] / r), min_height))
    r = height1 / min_height
    prefix_icon_image = prefix_icon_image.resize((int(prefix_icon_image.size[0] / r), min_height))


    (width1, height1) = prefix_icon_image.size
    (width2, height2) = image.size

    result_width = width1 + width2
    result_height = min(height1, height2)

    result = Image.new('RGBA', (result_width, result_height))
    result.paste(im=prefix_icon_image, box=(0, 0))
    result.paste(im=image, box=(width1, 0))

    return QtGui.QIcon(ImageQt.toqpixmap(result))



def convert_to_QPixmap_with_prefix_icon(text, prefix_icon=None, fs=9, darkmode=False, scale=2):
    qpixmap = convert_to_QPixMap(text, fs=fs, darkmode=darkmode, scale=scale)
    qicon = QtGui.QIcon(qpixmap)

    if prefix_icon is None:
        return qicon

    prefix_icon_image = ImageQt.fromqimage(prefix_icon.pixmap(prefix_icon.actualSize(QtCore.QSize(1024, 1024))).toImage())
    image = ImageQt.fromqimage(qicon.pixmap(qicon.actualSize(QtCore.QSize(1024, 1024))).toImage())


    (width1, height1) = prefix_icon_image.size
    (width2, height2) = image.size

    min_height = min(height1, height2)

    r = height2 / min_height
    image = image.resize((int(image.size[0] / r), min_height))
    r = height1 / min_height
    prefix_icon_image = prefix_icon_image.resize((int(prefix_icon_image.size[0] / r), min_height))


    (width1, height1) = prefix_icon_image.size
    (width2, height2) = image.size

    result_width = width1 + width2
    result_height = min(height1, height2)

    result = Image.new('RGBA', (result_width, result_height))
    result.paste(im=prefix_icon_image, box=(0, 0))
    result.paste(im=image, box=(width1, 0))
        
    result = ImageQt.toqpixmap(result)
    result.setDevicePixelRatio(scale)
    return result