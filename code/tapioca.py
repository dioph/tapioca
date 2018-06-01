import numpy as np
from astropy.io import ascii
import pyqtgraph as pg
from pyqtgraph import QtGui
from PyQt5.QtWidgets import *
import sys, re, urllib
from utils import *


def make_url(select='', where='', form='', order=''):
    url = "https://exoplanetarchive.ipac.caltech.edu/" \
          "cgi-bin/nstedAPI/nph-nstedAPI?table=q1_q17_dr25_stellar"
    columns = ['', 'st_delivname', 'kepid', 'tm_designation', 'ra', 'dec',
               'kepmag', 'teff', 'teff_err1', 'teff_err2', 'teff_prov',
               'logg', 'logg_err1', 'logg_err2', 'logg_prov', 'feh',
               'feh_err1', 'feh_err2', 'feh_prov', 'radius', 'radius_err1',
               'radius_err2', 'mass', 'mass_err1', 'mass_err2', 'dens',
               'dens_err1', 'dens_err2', 'prov_sec', 'nconfp', 'nkoi',
               'ntce', 'st_quarters', 'st_vet_date_str']
    formats = ['', 'ascii', 'ipac', 'bar', 'xml', 'json']

    w = re.split('<|=|>| and ', where)
    s = re.split(',', select)
    o = re.split(',| desc', order)
    for prop in s:
        assert prop in columns, "unknown property {}".format(prop)
    for prop in o:
        assert prop in columns, "unknown property {}".format(prop)
    for prop in w[::2]:
        assert prop in columns, "unknown property {}".format(prop)
    for val in w[1::2]:
        assert np.isscalar(val), "{} is not a scalar".format(val)

    assert form in formats, "unknown format {}".format(form)
    if select is not '':
        url += "&select={}".format(select)
    if where is not '':
        url += "&where={}".format(where)
    if form is not '':
        url += "&format={}".format(form)
    if order is not '':
        url += "&order={}".format(order)
    return url


class interface(QWidget):
    def __init__(self, parent=None):
        super(interface, self).__init__(parent=parent)
        self.lst_size = 10
        self.layout = QHBoxLayout(self)

        self.url = make_url()

        self.lay_search = QVBoxLayout()
        self.bt_search = QPushButton('SEARCH')
        self.bt_search.clicked.connect(self.search)
        self.lab_search_select = QLabel('SELECT')
        self.txt_search_select = QLineEdit()
        self.lab_search_where = QLabel('WHERE')
        self.txt_search_where = QLineEdit()
        self.lab_search_form = QLabel('FORMAT')
        self.txt_search_form = QLineEdit()
        self.lab_search_order = QLabel('ORDER')
        self.txt_search_order = QLineEdit()
        self.lay_search.addWidget(self.bt_search)
        self.lay_search.addWidget(self.lab_search_select)
        self.lay_search.addWidget(self.txt_search_select)
        self.lay_search.addWidget(self.lab_search_where)
        self.lay_search.addWidget(self.txt_search_where)
        self.lay_search.addWidget(self.lab_search_form)
        self.lay_search.addWidget(self.txt_search_form)
        self.lay_search.addWidget(self.lab_search_order)
        self.lay_search.addWidget(self.txt_search_order)
        self.msg_search = QLabel(self)
        self.lst_search = QListWidget()
        self.lst_search.itemActivated.connect(self.info)
        self.lay_search.addWidget(self.msg_search)
        self.lay_search.addWidget(self.lst_search)
        self.layout.addLayout(self.lay_search)

        self.lay_info = QVBoxLayout()
        self.lst_info = QListWidget()
        self.lay_info.addWidget(self.lst_info)
        self.layout.addLayout(self.lay_info)
        
        self.win = pg.GraphicsWindow()
        self.layout.addWidget(self.win)
        
        self.lay_menu = QVBoxLayout()
        self.menu = QComboBox()
        self.menu.addItem('')
        self.menu.addItem('TPF')
        self.menu.activated.connect(self.update)
        self.lay_menu.addWidget(self.menu)
        self.layout.addLayout(self.lay_menu)
        
    def search(self):
        select = self.txt_search_select.text()
        if 'kepid' not in select:
            select = 'kepid,'+select
        where = self.txt_search_where.text()
        form = self.txt_search_form.text()
        order = self.txt_search_order.text()
        self.lst_search.clear()
        self.msg_search.setText('')
        try:
            url = make_url(select, where, form, order)
            url = urllib.parse.quote(url, safe=':/&=?')
        except AssertionError as err:
            print(select, where, form)
            self.msg_search.setText('<font color="red">{}</font>'.format(str(err).upper()))
            return
        self.data = ascii.read(url)
        kics = self.data['kepid']
        for kic in kics[:self.lst_size]:
            self.lst_search.addItem(str(kic))
        
    def info(self, item):
        self.lst_info.clear()
        idx = self.lst_search.currentRow()
        star = self.data[idx]
        self.kic = star['kepid']
        for col in star.colnames:
            self.lst_info.addItem(u"{0}: {1}".format(col, star[col]))

    def update(self, arg_1):
        pass
           
def start():
    global app
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    font = QtGui.QFont('Courier')
    app.setFont(font)
    new = interface()
    new.show()
    app.exit(app.exec_())
    return new
    
if __name__ == "__main__":
    new = start()
