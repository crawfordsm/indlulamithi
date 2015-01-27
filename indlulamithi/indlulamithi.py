#!/usr/bin/env python

import sys
import os 
import shutil 
import time 
import datetime
import glob
import pickle
import logging
logging.disable('warning')
logging.disable('info')

import numpy as np
import matplotlib.pyplot as plt


from astropy.io import fits

from PyQt4 import QtGui, QtCore

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar


import  giraffe_reduce as gr
from makeorders import make_orders
from extract_spectra import find_first_order, fit_all_orders, extract_spectra, get_wavelength

def indlulamithi(datedir, obsdate=None, copy_data=False):
    """indlulamithi is the interactive data reduction program for the Giraffe
       spectragraph on the 1.9m telescope

    Parameters
    ----------
    datadir: str
        Directory containing the data to be reduced.
    """
    #create GUI
    App = QtGui.QApplication([])

    aw = indlulamithiWindow(datedir, obsdate=obsdate, copy_data=copy_data)

    aw.setMinimumHeight(900)
    aw.setMinimumWidth(900)
    aw.show()

    # Start application event loop
    exit=App.exec_()

    return

class indlulamithiWindow(QtGui.QMainWindow):

    def __init__(self, datadir,  obsdate=None, copy_data=False, run_clean=True):
        self.datadir = datadir
        self.obsdate = obsdate
        self.copy_data = copy_data
        if obsdate is not None:
            self.archivedir = '/data/74in/giraffe/data/image/%s/' % obsdate 
        self.w1=4500
        self.w2=5500

        #set up the observing log
        self.header_list = gr.default_headers + [' S/N ', ' n1   ', ' n2   ', ' w1      ', ' w2      ']
        print self.header_list
        self.image_list = glob.glob(datadir+'a*fits')
        self.obs_dict = gr.create_observing_dict(self.image_list, header_list=self.header_list)
        self.nrow = len(self.obs_dict)
        self.load_data()
  
      
        # Setup widget
        QtGui.QMainWindow.__init__(self)

        # Set main widget
        self.main = QtGui.QWidget(self)

        # Set window title
        self.setWindowTitle("Indlulamithi--Giraffe Data Reduction")

        #set up some buttons and labels
        self.calLabel = QtGui.QLabel("Calibrations:")
        self.calLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised )
        self.flatButton = QtGui.QPushButton("Make Camera Flat")
        self.flatButton.clicked.connect(self.run_make_camera_flat)
        self.proButton = QtGui.QPushButton("Process Data")
        self.proButton.clicked.connect(self.run_clean_data)
        self.arcButton = QtGui.QPushButton("Arc Drift")
        self.arcButton.clicked.connect(self.close)

        welcome_msg = logging.info('Welcome to Indlulamithi--Giraffe Data Reduction')
        self.msgBox = QtGui.QTextEdit(welcome_msg) 
        self.msgBox.setMaximumHeight(150)
        self.msgBox.setReadOnly(True)
      
        self.w1Label = QtGui.QLabel("W1:")
        self.w1Label.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised )
        self.w1ValueLabel = QtGui.QLineEdit(str(self.w1))
        self.w2Label = QtGui.QLabel("W2:")
        self.w2Label.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised )
        self.w2ValueLabel = QtGui.QLineEdit(str(self.w2))
        self.plotButton = QtGui.QPushButton('plot')
        self.plotButton.clicked.connect(self.plotspectra)
        
        self.quitButton = QtGui.QPushButton("Quit")
        self.quitButton.clicked.connect(self.close)



        self.obstable=QtGui.QTableWidget()
        self.obstable.setRowCount(self.nrow)
        self.obstable.setColumnCount(len(self.header_list))
        self.obstable.setHorizontalHeaderLabels(self.header_list)

        #if run_clean: self.run_clean_data()
        self.update_table()


        #set up the spectral plot
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)


        mainLayout = QtGui.QGridLayout(self.main)
        #add buttons for calibrating the data
        #mainLayout.addWidget(self.msgBox,0,0,2,5)

        mainLayout.addWidget(self.calLabel,1,0,1,1)
        mainLayout.addWidget(self.flatButton,1,1,1,1)
        mainLayout.addWidget(self.proButton,1,2,1,1)
        #mainLayout.addWidget(self.arcButton,1,3,1,1)

        mainLayout.addWidget(self.obstable,2,0,2,5)
        mainLayout.addWidget(self.canvas,4,0,2,5)
        mainLayout.addWidget(self.toolbar,6,0,1,5)

        mainLayout.addWidget(self.w1Label,7,0,1,1)
        mainLayout.addWidget(self.w1ValueLabel,7,1,1,1)
        mainLayout.addWidget(self.w2Label,7,2,1,1)
        mainLayout.addWidget(self.w2ValueLabel,7,3,1,1)
        mainLayout.addWidget(self.plotButton,7,4,1,1)

        mainLayout.addWidget(self.quitButton,8,0,1,5)
        
        self.setCentralWidget(self.main)

        # Destroy widget on close
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)

    def update_table(self):
        keys = self.obs_dict.keys()
        keys.sort()
        self.nrows = len(keys)
        print self.nrows, self.obstable.rowCount()
        for i in range(self.obstable.rowCount(), self.nrows):
            print i
            self.obstable.insertRow(i-1)
       
        for i, k in enumerate(keys):
           self.setrow(i, resize=True)


        for j in range(0, len(self.header_list)):
           self.obstable.resizeColumnToContents(j)

    def plotspectra(self):
        i = self.obstable.currentRow()
        specfile='s' + str(self.obstable.item(i,0).text())
        print specfile
        if not os.path.isfile(specfile):
           msg = 'Invalid row selected'
           print msg
           return 
        hdu = fits.open(specfile.strip())
        warr = hdu[1].data['wavelength']
        farr = hdu[1].data['flux']
        logging.info('Plotting Spectra for %s' % specfile)
        ax = self.figure.add_subplot(111)
        ax.hold(False)
        ax.plot(warr, farr)
        w1 = float(self.w1ValueLabel.text())
        w2 = float(self.w2ValueLabel.text())
        ax.axis([w1,w2, farr.min(), farr.max()])
        ax.set_xlabel('Wavelength ($\AA$)')
        ax.set_ylabel('Counts')
        self.canvas.draw()

    def run_clean_data(self):
        """Process all files"""

        if self.copy_data:
           print 'Running copy'
           img_list = glob.glob(self.archivedir+"a*fits")
           for img in img_list:
                 if not os.path.isfile(self.datadir+os.path.basename(img)):
                      shutil.copy(img, self.datadir)
             

        self.image_list = glob.glob(self.datadir+'a*fits')
        self.obs_dict = gr.create_observing_dict(self.image_list, header_list=self.header_list)
  
        #check to see if the camera flats exist
        if not os.path.isfile(self.datadir+'camera_flat.fits'):
            logging.warning('\ncamera_flat.fits does not exist.  Not processing data.')
            return

        keys = self.obs_dict.keys()
        keys.sort()
        #first run through flats
        for k in keys:
            if self.obs_dict[k][1]=='FLAT':
                 self.clean_data(k, self.obs_dict[k])
        #then the arcs
        for k in keys:
            if self.obs_dict[k][1]=='ARC':
                 self.clean_data(k, self.obs_dict[k])
        #then everything else
        for k in keys:
            if self.obs_dict[k][1] not in ['CAMERA', 'JUNK', 'ARC', 'FLAT']:
                 self.clean_data(k, self.obs_dict[k])
        self.load_data()
        self.update_table()
        print 'Finished Processing Data'

    def load_data(self):
        keys = self.obs_dict.keys()
        keys.sort()
        for img in keys:
            profile=self.datadir+'r'+img  
            if self.obs_dict[img][1]=='ARC':
                 solfile=self.datadir+'r'+img.replace('fits', 'pkl')
            elif self.obs_dict[img][1] in ['FLAT', 'CAMERA', 'JUNK']:
                 continue
            else:
                 img_arc = self.find_closest_image(img, 'ARC')
                 solfile = 'r' + img_arc.replace('fits', 'pkl')
                 try:
                     hdu = fits.open('s'+img)
                     try:
                          self.obs_dict[img][7] = '%i' % hdu[1].header['SN']
                     except:
                          pass
                     hdu.close()
                 except:
                     pass
                 
            if os.path.isfile(solfile):
                sol_dict = pickle.load(open(solfile))
                keys = np.array(sol_dict.keys())
                n1 = keys.min()
                n2 = keys.max()
                ws = sol_dict[n1][4]
                w2 = ws(24)
                ws = sol_dict[n2][4]
                w1 = ws(1077)
                self.obs_dict[img][8]='%3i' % int(n1)
                self.obs_dict[img][9]='%3i' % int(n2)
                self.obs_dict[img][10]='%6.2f' % w1
                self.obs_dict[img][11]='%6.2f' % w2


    def clean_data(self, img, hdr):
        """Fully process a single file"""
   
        #check to see if the camera flats exist
        cam_flat = self.datadir+'camera_flat.fits'
        if not os.path.isfile(cam_flat):
            logging.warning('\ncamera_flat.fits does not exist.  Not processing data.' % img)
            return

        profile=self.datadir+'r'+img    
        #first check to see if the data has already been reduced
        if os.path.isfile(profile): 
            logging.warning('Image %s has already been processed' % img)
            return

        #pass it through basic processing
        print 'Processing %s' % img
        logging.info('Processing %s' % img)
        ccd = gr.giraffe_reduce(self.datadir+img, camera_flat=cam_flat, outfile=profile)
 
        #check what type of data it is and then run any more
        #advance processing on it
        if hdr[1]=='FLAT': 
            ordfile=profile.replace('fits', 'orders')
            make_orders(ccd.data, xc=680, limit=1.5, image_size=10, order_size=2, outfile=ordfile)
        elif hdr[1]=='CAMERA' or hdr[1]=='JUNK':
            pass
        elif hdr[1]=='ARC':
            img_orders = self.find_closest_image(img, 'FLAT')
            if img_orders is None: 
                 logging.info('\nNo fiber flats are available.  Not processing data.' % img)
                 return 
            img_orders = 'r' + img_orders.replace('fits', 'orders')
            orders = np.loadtxt(img_orders)
            arcfile = os.path.dirname(gr.__file__)+'/data/thar.fits'
            n1 = find_first_order(ccd.data, orders, arcfile=arcfile, n1=55, n2=160)
            print 'First order: ', n1 
            print 'Last order: ', n1 + len(orders)
            print 'First Wavelength:', get_wavelength(680, n1+len(orders))
            print 'Last Wavelength:', get_wavelength(680, n1)
            solfile = profile.replace('fits', 'pkl')
            wsarcfile = os.path.dirname(gr.__file__)+'/data/thar_list.txt'
            fit_all_orders(ccd.data, orders, n1, outfile=solfile, 
                           arcfile=arcfile,  wsarcfile=wsarcfile, dw=3, 
                           gamma = 6.4, cam_foc=475)
            outfile = 's' + img
            sol_dict = pickle.load(open(solfile))
            xarr, warr, farr, farr_err, norders, sn = \
                  extract_spectra(ccd, sol_dict, dy=5, outfile=outfile)
        else:
            img_arc = self.find_closest_image(img, 'ARC')
            img_arc = 'r' + img_arc.replace('fits', 'pkl')
            sol_dict = pickle.load(open(img_arc))
            outfile = 's' + img
            xarr, warr, farr, farr_err, norders, sn = \
                  extract_spectra(ccd, sol_dict, dy=5, outfile=outfile)
 
        return 
        

    def find_closest_image(self, img, imtype):
        """For a given imtype, find the closest image in time to that
           image type
        """ 
        obstime = datetime.datetime.strptime('%s-%s-%s %s' % ('2015', '01', '23', self.obs_dict[img][3]), '%Y-%m-%d %H:%M:%S')
        if obstime.hour < 12: ot = obstime + datetime.timedelta(days=1)
        keys = self.obs_dict.keys()
        keys.sort()
        best_time=1e5
        best_image = None
        for k in keys:
            if self.obs_dict[k][1]==imtype:
                ot = datetime.datetime.strptime('%s-%s-%s %s' % ('2015', '01', '23', self.obs_dict[k][3]), '%Y-%m-%d %H:%M:%S')
                if ot.hour < 12: ot = ot + datetime.timedelta(days=1)
                t = min((obstime-ot).seconds, (ot-obstime).seconds)
                if t < best_time:
                   best_time =  t
                   best_img = k
        return best_img
        
        
    def run_make_camera_flat(self):
        """Create the camera flat"""
        outfile = 'camera_flat.fits'
        img_list = []
        for k in self.obs_dict:
            if self.obs_dict[k][1] == 'CAMERA': img_list.append(self.obs_dict[k][0])
        print img_list
        logging.info('\nCreate a new camera flat with name %s\n' % outfile) 
        gr.make_camera_flat(img_list, outfile)
        

    def setrow(self, i, resize=True):
        """Set all the values in a row from the obsdictionary"""
        if i >= len(self.obs_dict): return

	keys=self.obs_dict.keys()
        keys.sort()
        k = keys[i]

        nameItem=QtGui.QTableWidgetItem(k)
        self.obstable.setItem(i, 0, nameItem)
        self.obstable.resizeColumnToContents(0)
        for j in range(1, len(self.obs_dict[k])):
           item=self.parseItem(self.obs_dict[k][j])
           self.obstable.setItem(i, j, item)
           if resize: self.obstable.resizeColumnToContents(j)

    def parseItem(self, x):
       """Parse an object so it can be entered into the table"""
       if isinstance(x, str):
           return QtGui.QTableWidgetItem(x)
       elif isinstance(x, float):
           return QtGui.QTableWidgetItem('%f' % x)
       elif isinstance(x, int):
           return QtGui.QTableWidgetItem('%i' % x)
       return QtGui.QTableWidgetItem('')



if __name__=='__main__':

    import getopt
    try:
       opts, args = getopt.getopt(sys.argv[1:],"h:d:",["help","date="])
    except getopt.GetoptError:
       sys.exit(2)
 
    obsdate=None
    for o, a in opts:
       if o in ("-h", "--help"):
           print indlulamithi.__doc__
           sys.exit()
       if o in ("-d", "--date"):
           if len(args)==0:
               obsdate = str(a)
           if len(args)==1:
               obsdate = str(args[0])
           if len(args)>1:
               print indlulamithi.__doc__
               sys.exit(2)

    if obsdate is None : 
       datadir=os.getcwd()+'/'
       copy_data=False
    else:
       if not os.path.isdir(obsdate):
            os.mkdir(obsdate)
       os.chdir(obsdate)
       
       datadir=os.getcwd()+'/'
       copy_data=True
       
    print datadir
    indlulamithi(datadir, obsdate=obsdate, copy_data=copy_data)
   
