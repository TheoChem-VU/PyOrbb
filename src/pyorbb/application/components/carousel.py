from PySide6 import QtWidgets, QtGui, QtCore
import webbrowser
import requests
from bs4 import BeautifulSoup
import tcmu
import datetime
import json
import pprint
from pyorbb.application.components import shadow


class FadeWidget(QtWidgets.QWidget):
    def __init__(self, side, parent=None):
        super().__init__(parent)
        self.side = side
        self.setAttribute(QtCore.Qt.WA_TransparentForMouseEvents)
        self.setAttribute(QtCore.Qt.WA_NoSystemBackground)
        self.update_background_color()

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        if self.side == 'right':
            gradient = QtGui.QLinearGradient(self.width(), 0, 0, 0)
        elif self.side == 'left':
            gradient = QtGui.QLinearGradient(0, 0, self.width(), 0)
        elif self.side == 'top':
            gradient = QtGui.QLinearGradient(0, self.height(), 0, 0)
        elif self.side == 'bottom':
            gradient = QtGui.QLinearGradient(0, 0, 0, self.height())

        solid = QtGui.QColor(self.background_color)
        solid.setAlpha(255)

        transparent = QtGui.QColor(self.background_color)
        transparent.setAlpha(0)

        gradient.setColorAt(0, solid)
        gradient.setColorAt(1, transparent)

        painter.fillRect(self.rect(), gradient)

    def update_background_color(self):
        darkmode = QtWidgets.QApplication.instance().isDarkMode
        self.background_color = '#3E3E3E' if darkmode else '#F7F7F7'


@tcmu.cache_file("cited_by.json", datetime.timedelta(weeks=1))
def _get_citedby_data(url: str):
    headers = {
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
        "AppleWebKit/537.36 (KHTML, like Gecko) "
        "Chrome/126.0 Safari/537.36"
        ),
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    "Accept-Language": "en-US,en;q=0.9",
        }
    res = requests.get(url, headers=headers, allow_redirects=True, timeout=1.5).content
    # with open('/Users/yumanhordijk/Desktop/test.html') as inp:
    #     res = inp.read()
    soup = BeautifulSoup(res, 'html.parser')
    # print(soup)
    # current_art = soup.find(id='gs_res_ccl_top')
    # current_art_title = current_art.find('h2').a.get_text()
    # current_art_title_doi = tcmu._get_doi_data_from_title(current_art_title)['DOI']
    # print(soup)
    citedby_data = []
    matches = soup.find_all('div', 'gs_ri')
    # print(matches)
    # print(matches)
    for match in matches:
        gs_a = match.find('div', 'gs_a')
        text = gs_a.get_text().replace('…', '')
        parts = text.split('- ')
        if len(parts) == 3:
            authors, journal_year, publisher = parts
            journal = journal_year.split(',')[0].strip()
            if journal.isnumeric():
                journal = None
        elif len(parts) == 2:
            authors, publisher = parts
            journal = None

        if journal is not None and 'chemrxiv' in journal:
            continue

        authors = authors.strip().split(',')
        authors = [author.strip() for author in authors]
        datum = {
            'link': match.h3.a.get('href'), 
            'title': match.h3.a.get_text(),
            'authors': authors,
            # 'journal': journal,
            }
        citedby_data.append(datum)

    citation_data = []
    for citedby in citedby_data:
        # print(citedby['title'])

        citation_datum = tcmu._get_doi_data_from_query(
            title=citedby['title'],
            author=citedby['authors'],)
            # container_title=citedby['journal'])

        # print(citation_datum['title'])
        # print()
        if citation_datum is None:
            continue
        # print(citation_datum)
        if 'reference' not in citation_datum:
            continue

        references = citation_datum["reference"]
        for reference in references:
            if 'DOI' not in reference:
                continue
            # if reference['DOI'] != current_art_title_doi:
            #     continue
            citation_data.append(citation_datum)
            break

    return citation_data



# # url = "http://webcache.googleusercontent.com/search?q=cache:https://scholar.google.com/scholar?hl=en&num=20&as_sdt=2005&sciodt=0,5&cites=13608188881403064766&scipsc=&q=&scisbd=1"
# url = "https://www.scopus.com/pages/publications/85208784619#tab=citedBy"
# print(_get_citedby_data(url))
# # print(_get_citedby_data(url))
# # print(_get_citedby_data(url))
# # print(_get_citedby_data(url))
# exit()

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

        self.setObjectName("Pub_main")
        self.setup()
        shadow.apply(self, radius=10)

        # grad = QtGui.QGradient(QtGui.QGradient.Preset.WarmFlame)
        # print(grad)
        # gradient = QLinearGradient(0, 0, 0, widget.height())  # top to bottom
        # gradient.setColorAt(0.0, QColor("#4facfe"))  # start color
        # gradient.setColorAt(1.0, QColor("#00f2fe"))  # end color

        # palette = self.palette()
        # palette.setBrush(QtGui.QPalette.ColorRole.Window, QtGui.QBrush(grad))
        # self.setPalette(palette)
        # self.setAutoFillBackground(True)

        self.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)

        self.setToolTip(self._link)

    # def paintEvent(self, event):
    #     super().paintEvent(event)
    #     painter = QtGui.QPainter(self)

    #     preset = QtGui.QGradient(QtGui.QGradient.Preset.ShadyWater)
    #     stops = preset.stops()

    #     gradient = QtGui.QLinearGradient(0, 0, 0, 1)  # top to bottom
    #     gradient.setCoordinateMode(QtGui.QGradient.CoordinateMode.ObjectMode)
    #     gradient.setStops(stops)

    #     painter.fillRect(self.rect(), gradient)


    def mousePressEvent(self, event):
        import webbrowser
        ret = webbrowser.open(self._link)

    def setup(self):
        month_name = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun', 7: 'Jul', 8: 'Aug', 9: 'Sept', 10: 'Oct', 11: 'Nov', 12: 'Dec'}
        journal = tcmu._get_journal_abbreviation(self._journal)
        title_label = QtWidgets.QLabel('<b>' + self._title + '</b>')

        title_label.setWordWrap(True)

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
        names = [f"<i>{last}</i>" for first, last in zip(initials, last_names)]
        if len(names) > 4:
            author_str = ', '.join(names[:3]) + ', <i>et al.</i>'
        else:
            author_str = ', '.join(names)

        author_label = QtWidgets.QLabel(author_str)
        author_label.setWordWrap(True)

        ref_label = QtWidgets.QLabel(f'<i><b>{journal}</b></i><br>{month_name[self._date[1]]} {self._date[0]}')

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
            title_label.setObjectName('carousel_title')
            _layout.addWidget(title_label)


        self.fade_left = FadeWidget(side='left', parent=self.scroll_area)
        self.fade_right = FadeWidget(side='right', parent=self.scroll_area)
        # self.fade_top = FadeWidget(side='top', parent=self.scroll_area)
        # self.fade_bottom = FadeWidget(side='bottom', parent=self.scroll_area)

        self.fade_width = 35

        # Connect scroll updates
        scrollbar = self.scroll_area.horizontalScrollBar()

        _layout.addWidget(self.scroll_area)

        self.setup()
        self.resize_fades()

    def update_background_color(self):
        self.fade_left.update_background_color()
        self.fade_right.update_background_color()
        # self.fade_top.update_background_color()
        # self.fade_bottom.update_background_color()
        ...

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.resize_fades()

    def resize_fades(self):
        # Make scroll area fill entire container
        # Position fade overlays ON TOP
        self.fade_left.setGeometry(
            120, 
            0, 
            self.fade_width, 
            self.scroll_area.height()
            )
        self.fade_right.setGeometry(
            self.scroll_area.width() - self.fade_width - 120,
            0,
            self.fade_width,
            self.scroll_area.height()
            )
        # self.fade_top.setGeometry(
        #     120, 
        #     self.scroll_area.height() - self.fade_width, 
        #     self.scroll_area.width() - 240, 
        #     self.fade_width
        #     )
        # self.fade_bottom.setGeometry(
        #     120, 
        #     0,
        #     self.scroll_area.width() - 240, 
        #     self.fade_width
        #     )

        # Ensure fades are above scroll area
        self.fade_left.raise_()
        self.fade_right.raise_()
        # self.fade_top.raise_()
        # self.fade_bottom.raise_()
        ...

    def setup(self):
        ...


class PublicationCarousel(Carousel):
    def __init__(self, parent):
        super().__init__(parent=parent, title='<b>Recent Publications Citing PyOrbb</b>')

    def setup(self):
        url = "https://scholar.google.com/scholar?cites=8000893946037734095&scisbd=1"
        # try to get data
        try:
            data = _get_citedby_data(url)
        # otherwise we display a 404 error message
        except:
            self._carousel_layout.addWidget(QtWidgets.QLabel('Sorry! Could not find the right data.'))
            return

        for row in data:
            try:
                self._carousel_layout.addWidget(PublicationWidget(self, row))
            except:
                pass

        self._carousel_layout.addStretch()
