import geemap
import ee
import math

class slope_correction(object):
    def __init__(self, collection):
        self.model = 'volume'
        self.elevation = ee.Image('USGS/SRTMGL1_003')
        self.buffer = 10
        self.ninetyRad = ee.Image.constant(90).multiply(math.pi/180)
        self.collection=collection

    def _volume_model(self,theta_iRad, alpha_rRad):
        nominator = (self.ninetyRad.subtract(theta_iRad).add(alpha_rRad)).tan()
        denominator = (self.ninetyRad.subtract(theta_iRad)).tan()
        return nominator.divide(denominator)

    def _surface_model(self,theta_iRad, alpha_rRad, alpha_azRad):
        nominator = (self.ninetyRad.subtract(theta_iRad)).cos()
        denominator = alpha_azRad.cos().multiply((self.ninetyRad.subtract(theta_iRad).add(alpha_rRad)).cos())
        return nominator.divide(denominator)

    def _erode(self,img, distance):
        d = (img.Not().unmask(1).fastDistanceTransform(30).sqrt().multiply(ee.Image.pixelArea().sqrt()))
        return img.updateMask(d.gt(distance))

    def _masking(self,alpha_rRad, theta_iRad, proj, buffer):
        layover = alpha_rRad.lt(theta_iRad).rename('layover')
        shadow = alpha_rRad.gt(ee.Image.constant(-1).multiply(self.ninetyRad.subtract(theta_iRad))).rename('shadow')
        mask = layover and shadow
        if (buffer > 0):
            mask = self._erode(mask, buffer)
        return mask.rename('no_data_mask')

    def _correct(self,image):
        geom = image.geometry()
        proj = image.select(1).projection()
        heading = (ee.Terrain.aspect(image.select('angle')).reduceRegion(ee.Reducer.mean(), geom, 1000).get('aspect'))
        sigma0Pow = ee.Image.constant(10).pow(image.divide(10.0))
        theta_iRad = image.select('angle').multiply(math.pi/180).clip(geom)
        phi_iRad = ee.Image.constant(heading).multiply(math.pi/180)
        alpha_sRad = ee.Terrain.slope(self.elevation).select('slope').multiply(math.pi/180).setDefaultProjection(proj).clip(geom)
        phi_sRad = ee.Terrain.aspect(self.elevation).select('aspect').multiply(math.pi/180).setDefaultProjection(proj).clip(geom)
        phi_rRad = phi_iRad.subtract(phi_sRad)
        alpha_rRad = (alpha_sRad.tan().multiply(phi_rRad.cos())).atan()
        alpha_azRad = (alpha_sRad.tan().multiply(phi_rRad.sin())).atan()
        gamma0 = sigma0Pow .divide(theta_iRad.cos())
        if (self.model == 'volume'):
          corrModel = self._volume_model(theta_iRad, alpha_rRad)
        if (self.model == 'surface'):
          corrModel = self._surface_model(theta_iRad, alpha_rRad, alpha_azRad)
        gamma0_flat = gamma0.divide(corrModel);
        gamma0_flatDB = (ee.Image.constant(10).multiply(gamma0_flat.log10()).select(['VV', 'VH']))
        mask = self._masking(alpha_rRad, theta_iRad, proj, self.buffer)
        return gamma0_flatDB.addBands(mask).copyProperties(image)
    def main(self):
        return self.collection.map(self._correct)