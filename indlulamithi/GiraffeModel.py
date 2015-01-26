import numpy as np
from PySpectrograph.Spectrograph import Spectrograph, Grating, Optics, CCD, Detector, Slit
from PySpectrograph import SpectrographError


class GiraffeError(Exception): 
    pass

class GiraffeModel (Spectrograph):
   """GiraffeModel is a class that describes the Giraffe Spectrograph located on the 1.9 m telescope

"""

   def __init__(self, grating_name='Giraffe', slit=2.0,  \
                      order=83, gamma=6.0, xbin=1, ybin=1, xpos=0.00, ypos=0.00):
       
       #set up the parts of the grating
       self.grating_name=grating_name


       #set the telescope
       self.set_telescope('1.9m')

       #set the collimator
       self.set_collimator('Giraffe')

       #set the camera
       self.set_camera('Giraffe')

       #set the detector
       self.set_detector('Giraffe', xbin=xbin, ybin=ybin, xpos=xpos, ypos=ypos)

       #set up the grating
       self.set_grating(self.grating_name, order=order)
 
       #set up the slit
       self.set_slit(slit)

       #set up the grating angle
       self.gamma=gamma

   def alpha(self, da=0.00):
       """Return the value of alpha for the spectrograph"""
       return self.grating.blaze+self.gamma 

   def beta(self, db=0):
       """Return the value of beta for the spectrograph
          
          Beta_o=(1+fA)*(camang)-gratang+beta_o
       """
       return self.grating.blaze-self.gamma+db

   def get_wavelength(self, xarr, gamma=0.0):
       """For a given spectrograph configuration, return the wavelength coordinate
          associated with a pixel coordinate.   

          xarr: 1-D Array of pixel coordinates
          gamma: Value of gamma for the row being analyzed

          returns an array of wavelengths in mm
       """
       d=self.detector.xbin*self.detector.pix_size*(xarr-self.detector.get_xpixcenter())
       dbeta=np.degrees(np.arctan(d/self.camera.focallength))
       return self.calc_wavelength(self.alpha(), -self.beta()+dbeta, gamma=gamma)
      


   def set_telescope(self, name='1.9m'):
       if name=='1.9m':
           self.telescope=Optics(name=name, focallength=34200.0)
       else:
           raise SpectrographError('%s is not a supported Telescope' % name)

   def set_collimator(self, name='Giraffe', focallength=447.0):
       if name=='Giraffe':
           self.collimator=Optics(name=name, focallength=focallength)
       else:
           raise SpectrographError('%s is not a supported collimator' % name)

   def set_camera(self, name='Giraffe', focallength=None):
       if name=='Giraffe':
           self.camera=Optics(name=name, focallength=400)
       else:
           raise SpectrographError('%s is not a supported camera' % name)

   


   def set_detector(self, name='Giraffe', geom=None, xbin=1, ybin=1, xpos=0, ypos=0):
       if name=='Giraffe':
               ccd=CCD(name='CCD1',  xpix=1128, ypix=1024, pix_size=0.024, xpos=0.00, ypos=0.00)
               self.detector=Detector(name=name, ccd=[ccd], xbin=xbin, ybin=ybin, \
                                      xpos=xpos, ypos=ypos)
       else:
           raise SpectrographError('%s is not a supported detector' % name)

   def set_grating(self, name=None, order=83):
       if name=='Giraffe':
          self.grating=Grating(name='Giraffe', spacing=31.6, blaze=65.5, order=order)
	  self.set_order(order)
       else:
          raise SpectrographError('%s is not a supported grating' % name)
      
   def set_order(self, order):
       self.order=order
       self.grating.order = order

   def set_slit(self, slitang=2.0):
       self.slit=Slit(name='Fiber', phi=slitang)
       self.slit.width=self.slit.calc_width(self.telescope.focallength)
