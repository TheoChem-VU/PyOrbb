from PySide6 import QtWidgets, QtGui, QtCore
import webbrowser
import requests
from bs4 import BeautifulSoup
import tcutility
import datetime
import json
import pprint
from pyfmo.application.components import shadow


class FadeWidget(QtWidgets.QWidget):
    def __init__(self, left_side=True, parent=None):
        super().__init__(parent)
        self.left_side = left_side
        self.setAttribute(QtCore.Qt.WA_TransparentForMouseEvents)
        self.setAttribute(QtCore.Qt.WA_NoSystemBackground)

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.Antialiasing)

        gradient = QtGui.QLinearGradient()

        if self.left_side:
            gradient = QtGui.QLinearGradient(self.width(), 0, 0, 0)
        else:
            gradient = QtGui.QLinearGradient(0, 0, self.width(), 0)

        solid = QtGui.QColor('#F8F8F8')
        solid.setAlpha(255)

        transparent = QtGui.QColor('#F8F8F8')
        transparent.setAlpha(0)

        gradient.setColorAt(0, solid)
        gradient.setColorAt(1, transparent)

        painter.fillRect(self.rect(), gradient)



@tcutility.cache_file("cited_by.json", datetime.timedelta(weeks=1))
def _get_citedby_data(url: str):
    headers = {"User-Agent": "Mozilla/5.0"}
    res = requests.get(url, headers=headers, allow_redirects=True).content
    soup = BeautifulSoup(res, 'html.parser')

    current_art = soup.find(id='gs_res_ccl_top')
    current_art_title = current_art.find('h2').a.get_text()
    current_art_title_doi = tcutility._get_doi_data_from_title(current_art_title)['DOI']

    match = soup.find(id='gs_res_ccl_mid')

    citedby_data = []
    for h3 in match.find_all('h3'):
        a = h3.a
        datum = {'link': a.get('href'), 'title': a.get_text()}
        citedby_data.append(datum)

    citation_data = []
    for citedby in citedby_data:
        citation_datum = tcutility._get_doi_data_from_title(citedby["title"])
        if 'reference' not in citation_datum:
            continue

        references = citation_datum["reference"]
        for reference in references:
            if 'DOI' not in reference:
                continue
            if reference['DOI'] != current_art_title_doi:
                continue
            citation_data.append(citation_datum)
            break

    return citation_data



class PublicationWidget(QtWidgets.QFrame):
    def __init__(self, parent, data):
        super().__init__(parent=parent)
        self.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.parent = parent
        self.layout = QtWidgets.QVBoxLayout(self)
        self._link = data["resource"]["primary"]["URL"]
        self._journal = data["container-title"][0]
        self._title = "".join("".join(data["title"]).splitlines())
        self._authors = data["author"]
        self._date = data["created"]["date-parts"][0]

        self.setObjectName("main")
        self.setStyleSheet("""
            QFrame#main {
                    border: none;
                    margin: 5px;
                    padding: 10px;
                    border-radius: 10px;
                    background-color: white;
                }
            QFrame#main:hover {
                background-color: #f6f6f6;
            }
            QLabel#Pub_ref {
                    font: 12px;
                    color: #1f2022;
                }
            QLabel#Pub_authors {
                    font: 12px;
                    color: #1f2022;
                }
            QLabel#Pub_title {
                    font: 12px;
                    color: #1f2022;
                }
            """)
        self.setup()
        # self.setMinimumWidth(250)
        shadow.apply(self)
        self.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)

    def mousePressEvent(self, event):
        import webbrowser
        # print(webbrowser.get('firefox'))
        ret = webbrowser.get('firefox').open(self._link)
        # ret = webbrowser.open(self._link)
        # print(ret)

    def setup(self):
        month_name = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun', 7: 'Jul', 8: 'Aug', 9: 'Sept', 10: 'Oct', 11: 'Nov', 12: 'Dec'}
        journal = tcutility._get_journal_abbreviation(self._journal)
        title_label = QtWidgets.QLabel('<b>' + self._title + '</b>')

        title_label.setWordWrap(True)
        title_label.setObjectName('Pub_title')

        initials = []
        last_names = []
        for author in self._authors:
            # we get the capital letters from the first names
            # these will become the initials for this author
            firsts = [char + "." for char in author["given"].title() if char.isupper()]
            firsts = " ".join(firsts)
            initials.append(firsts)
            last_names.append(author["family"].title())

        # format the citation correctly
        names = [f"{first} {last}" for first, last in zip(initials, last_names)]
        if len(names) > 3:
            author_str = ',  '.join(names[:3]) + ', <i>et al.</i>'
        else:
            author_str = ',  '.join(names)

        author_label = QtWidgets.QLabel(author_str)
        author_label.setWordWrap(True)
        author_label.setObjectName('Pub_authors')

        ref_label = QtWidgets.QLabel(f'<i><b>{journal}</b></i><br>{month_name[self._date[1]]} {self._date[0]}')
        # ref_label.setWordWrap(True)
        ref_label.setObjectName('Pub_ref')

        self.layout.addWidget(ref_label)
        self.layout.addWidget(title_label)
        self.layout.addWidget(author_label)
        self.layout.addStretch()


class Carousel(QtWidgets.QWidget):
    def __init__(self, parent, title=None):
        super().__init__(parent=parent)
        self.parent = parent
        self._items = []

        _layout = QtWidgets.QVBoxLayout(self)

        self._next_btn = QtWidgets.QPushButton('>')
        self._carousel_layout = QtWidgets.QHBoxLayout()

        self.scroll_area = QtWidgets.QScrollArea(self)
        self.scroll_area.setWidgetResizable(True)
        self.scroll_area_contents_widget = QtWidgets.QWidget()
        self.scroll_area.setWidget(self.scroll_area_contents_widget)
        self.scroll_area_contents_widget.setLayout(self._carousel_layout)
        self._prev_btn = QtWidgets.QPushButton('<')

        self.scroll_area.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.scroll_area.verticalScrollBar().setEnabled(False)
        self.scroll_area.horizontalScrollBar().setEnabled(True)
        self.scroll_area.setObjectName('carousel')

        if title is not None:
            title_label = QtWidgets.QLabel(title)
            title_label.setObjectName('title')
            _layout.addWidget(title_label)


        self.fade_left = FadeWidget(left_side=False, parent=self.scroll_area)
        self.fade_right = FadeWidget(left_side=True, parent=self.scroll_area)

        self.fade_width = 20

        # Connect scroll updates
        scrollbar = self.scroll_area.horizontalScrollBar()

        _layout.addWidget(self.scroll_area)

        if self.parent.parent.isDarkMode:
            self.setStyleSheet('''
                QLabel#title {
                    font-size: 20px;
                    color: #b0b0b0;
                }
                QWidget {
                    border: none;
                    background-color: transparent;
                    border-radius: 5px;
                }
                QPushButton {
                    font-size: 12px;
                    border: 1px solid darkgray;
                    border-radius: 13px;
                    padding: 8px;
                    margin: 0px; 
                    background-color: #3A3A3A;
                }
                QPushButton:hover {
                    background-color: #656565;
                }
                 QScrollBar:horizontal
                {
                    height: 15px;
                    margin: 3px 15px 3px 15px;
                    border: 1px transparent #2A2929;
                    border-radius: 2px;
                    background-color: white;    /* #2A2929; */
                }

                QScrollBar::handle:horizontal
                {
                    background-color: lightgrey;      /* #605F5F; */
                    min-width: 5px;
                    border-radius: 2px;
                }

                QScrollBar::add-line:horizontal
                {
                    margin: 0px 3px 0px 3px;
                    border-image: url(:/qss_icons/rc/right_arrow_disabled.png);
                    width: 10px;
                    height: 10px;
                    subcontrol-position: right;
                    subcontrol-origin: margin;
                }

                QScrollBar::sub-line:horizontal
                {
                    margin: 0px 3px 0px 3px;
                    border-image: url(:/qss_icons/rc/left_arrow_disabled.png);
                    height: 10px;
                    width: 10px;
                    subcontrol-position: left;
                    subcontrol-origin: margin;
                }

                QScrollBar::add-line:horizontal:hover,QScrollBar::add-line:horizontal:on
                {
                    border-image: url(:/qss_icons/rc/right_arrow.png);
                    height: 10px;
                    width: 10px;
                    subcontrol-position: right;
                    subcontrol-origin: margin;
                }


                QScrollBar::sub-line:horizontal:hover, QScrollBar::sub-line:horizontal:on
                {
                    border-image: url(:/qss_icons/rc/left_arrow.png);
                    height: 10px;
                    width: 10px;
                    subcontrol-position: left;
                    subcontrol-origin: margin;
                }

                QScrollBar::up-arrow:horizontal, QScrollBar::down-arrow:horizontal
                {
                    background: none;
                }


                QScrollBar::add-page:horizontal, QScrollBar::sub-page:horizontal
                {
                    background: none;
                }

            ''')
        else:
            self.setStyleSheet('''  
                QLabel#title {
                    font-size: 20px;
                    color: #808080;
                    padding-top: 10px;
                    padding-left: 90px;
                }  
                QLabel#ref {
                    font: 20px black;
                }  
                QLabel#authors {
                    font: 20px black;
                }  
                QLabel#article_title {
                    font: 20px black;
                }
                QWidget {
                    border: none;
                    background-color: transparent;
                    border-radius: 5px;
                }  
                QWidget#carousel {
                    padding-left: 0px;
                    padding-right: 0px;
                }  
                QPushButton {
                    font-size: 12px;
                    border: 1px solid lightgray;
                    border-radius: 13px;
                    padding: 8px;
                    margin: 0px; 
                    background-color: white;
                }
                QPushButton:hover {
                    background-color: #f0f0f0;
                }
                 QScrollBar:horizontal
                {
                    height: 11px;
                    margin: 3px 15px 3px 15px;
                    border: 1px transparent #2A2929;
                    border-radius: 4px;
                    background-color: transparent;    /* #2A2929; */
                }

                QScrollBar::handle:horizontal
                {
                    background-color: lightgrey;      /* #605F5F; */
                    min-width: 5px;
                    border-radius: 4px;
                }

                QScrollBar::add-line:horizontal
                {
                    margin: 0px 3px 0px 3px;
                    border-image: url(:/qss_icons/rc/right_arrow_disabled.png);
                    width: 10px;
                    height: 10px;
                    subcontrol-position: right;
                    subcontrol-origin: margin;
                }

                QScrollBar::sub-line:horizontal
                {
                    margin: 0px 3px 0px 3px;
                    border-image: url(:/qss_icons/rc/left_arrow_disabled.png);
                    height: 10px;
                    width: 10px;
                    subcontrol-position: left;
                    subcontrol-origin: margin;
                }

                QScrollBar::add-line:horizontal:hover,QScrollBar::add-line:horizontal:on
                {
                    border-image: url(:/qss_icons/rc/right_arrow.png);
                    height: 10px;
                    width: 10px;
                    subcontrol-position: right;
                    subcontrol-origin: margin;
                }


                QScrollBar::sub-line:horizontal:hover, QScrollBar::sub-line:horizontal:on
                {
                    border-image: url(:/qss_icons/rc/left_arrow.png);
                    height: 10px;
                    width: 10px;
                    subcontrol-position: left;
                    subcontrol-origin: margin;
                }

                QScrollBar::up-arrow:horizontal, QScrollBar::down-arrow:horizontal
                {
                    background: none;
                }


                QScrollBar::add-page:horizontal, QScrollBar::sub-page:horizontal
                {
                    background: none;
                }

            ''')

        self.setup()
        self.resize_fades()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.resize_fades()

    def resize_fades(self):
        # Make scroll area fill entire container
        # Position fade overlays ON TOP
        self.fade_left.setGeometry(120, 0, self.fade_width, self.scroll_area.height())
        self.fade_right.setGeometry(
            self.scroll_area.width() - self.fade_width - 120,
            0,
            self.fade_width,
            self.scroll_area.height()
        )

        # Ensure fades are above scroll area
        self.fade_left.raise_()
        self.fade_right.raise_()

    def setup(self):
        ...


class PublicationCarousel(Carousel):
    def __init__(self, parent):
        super().__init__(parent=parent, title='<b>Recent Publications Citing PyOrbb</b>')
        # shadow.apply(self, radius=-40)

    def setup(self):
        url = "https://scholar.google.com/scholar?hl=nl&as_sdt=2005&sciodt=0,5&cites=8000893946037734095&scipsc=&q=&scisbd=1"
        # url = "bla"
        # try to get data
        try:
            data = _get_citedby_data(url)
        # otherwise we display a 404 error message
        except:
            self._carousel_layout.addWidget(QtWidgets.QLabel('Sorry! Could not find the right data.'))

            return

        for row in data:
            # print(row)
            try:
                self._carousel_layout.addWidget(PublicationWidget(self, row))
            except:
                pass

        self._carousel_layout.addStretch()
