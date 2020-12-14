import scipy
from scipy.stats import rv_continuous,gaussian_kde,norm as normal
import numpy as np
import os

__all__=['get_distribution']

class skewed_normal(rv_continuous):
    "Skewed Normal Distribution"
    def _pdf(self,x,mu,left_sigma,right_sigma):
        try:
                mu=list(mu)[0]
                left_sigma=list(left_sigma)[0]
                right_sigma=list(right_sigma)[0]
        except:
            pass

        left=normal(loc=mu,scale=left_sigma)
        right=normal(loc=mu,scale=right_sigma)
        pdf=np.piecewise(x,[x<mu,x>=mu],
                            [lambda y : left.pdf(y)/np.max(left.pdf(y)),
                     lambda y : right.pdf(y)/np.max(right.pdf(y))])
        return(pdf/np.sum(pdf))

    def _argcheck(self,*args):
        return True

def _skewed_normal(name,dist_dat,dist_type):
    if dist_type+'_DIST_LIMITS' in dist_dat:
        a,b=dist_dat[dist_type+'_DIST_LIMITS']
    else:
        a=dist_dat[dist_type+'_DIST_PEAK']-3*dist_dat[dist_type+'_DIST_SIGMA'][0]
        b=dist_dat[dist_type+'_DIST_PEAK']+3*dist_dat[dist_type+'_DIST_SIGMA'][1]
    if a==b:
            return(lambda :[a])
    dist = skewed_normal(name,a=a,b=b)
    sample=np.linspace(a,b,int(1e4))
    return(lambda : np.random.choice(sample,1,
                                   p=dist._pdf(sample,dist_dat[dist_type+'_DIST_PEAK'],
                                    dist_dat[dist_type+'_DIST_SIGMA'][0],dist_dat[dist_type+'_DIST_SIGMA'][1])))


def _param_from_dist(dist_file,path):
    dist=np.loadtxt(os.path.join(path,dist_file))
    a=np.min(dist)-abs(np.min(dist))
    b=np.max(dist)+abs(np.max(dist))
    sample=np.linspace(a,b,int(1e4))
    pdf=gaussian_kde(dist.T).pdf(np.linspace(a,b,int(1e4)))
    return(lambda : np.random.choice(sample,1,p=pdf/np.sum(pdf)))

def _get_zdepend(dist_file,path,typ):
    fv=0 if typ=='add' else 1
    z,factor=np.loadtxt(os.path.join(path,dist_file),unpack=True)
    return(interp1d(z,factor,fill_value=fv,bounds_error=False))

def get_distribution(name,dist_dat,path,sn_or_host):
    dist_dict={}
    param=True
    if np.any(['DIST_FILE' in x for x in dist_dat.keys() if 'PARAM' in x]):

        if sn_or_host+'_PARAM_DIST_FILE' in dist_dat.keys():
            param_func=_param_from_dist(dist_dat[sn_or_host+'_PARAM_DIST_FILE'],path)
        else:
                raise RuntimeError("You may have a typo, did you mean to set 'PARAM_DIST_FILE' for %s?"%sn_or_host)


    elif np.any(['PARAM' in x for x in dist_dat.keys()]):
        try:
            param_func=_skewed_normal(name,dist_dat,sn_or_host+'_PARAM')
        except:
            raise RuntimeError("You may have a typo in the variables of your %s param distribution(s)."%sn_or_host)
    else:
        param=False

    if np.any(['DIST_FILE' in x for x in dist_dat.keys() if 'SCALE' in x]):
        if 'SCALE_DIST_FILE' in dist_dat.keys():
            scale_func=_param_from_dist(dist_dat['SCALE_DIST_FILE'],path)
        else:
            raise RuntimeError("You may have a typo, did you mean to set 'SCALE_DIST_FILE'?")

    elif np.any(['SCALE' in x for x in dist_dat.keys()]):
        try:
            scale_func=_skewed_normal(name,dist_dat,'SCALE')
        except:
            raise RuntimeError("You may have a typo in the variables of your %s scale distribution(s)."%sn_or_host)
    else:
        raise RuntimeError("Must supply scale distribution for every effect.")

    if np.any(['ZDEPEND' in x for x in dist_dat.keys() if 'PARAM' in x]):
        if sn_or_host+'_PARAM_ZDEPEND_FILE' in dist_dat.keys():

            if sn_or_host+'_PARAM_ZDEPEND_TYPE' in dist_dat.keys() and dist_dat[sn_or_host+'_PARAM_ZDEPEND_TYPE']=='MULTIPLY':
                zdepend=_get_zdepend(dist_dat[sn_or_host+'_PARAM_ZDEPEND_FILE'],path,'multiply')
                dist_dict['PARAM']=lambda z: (param_func()[0])*zdepend(z)
            else:
                zdepend=_get_zdepend(dist_dat[sn_or_host+'_PARAM_ZDEPEND_FILE'],path,'add')
                dist_dict['PARAM']=lambda z: param_func()[0]+zdepend(z)
        else:
            raise RuntimeError("You may have a typo, did you mean to set 'ZDEPEND' for %s?"%sn_or_host)
    elif param:
        dist_dict['PARAM']=lambda z:param_func()[0]

    if np.any(['ZDEPEND' in x for x in dist_dat.keys() if 'SCALE' in x]):
        if sn_or_host+'_SCALE_ZDEPEND_FILE' in dist_dat.keys():
            if sn_or_host+'_SCALE_ZDEPEND_TYPE' in dist_dat.keys() and\
               dist_dat[sn_or_host+'_SCALE_ZDEPEND_TYPE'].upper()=='MULTIPLY':
                zdepend=_get_zdepend(dist_dat[sn_or_host+'_SCALE_ZDEPEND_FILE'],path,'multiply')
                dist_dict['SCALE']=lambda z: (scale_func()[0])*zdepend(z)
            else:
                zdepend=_get_zdepend(dist_dat[sn_or_host+'_SCALE_ZDEPEND_FILE'],path,'add')
                dist_dict['SCALE']=lambda z: scale_func()[0]+zdepend(z)
        else:
            raise RuntimeError("You may have a typo, did you mean to set 'ZDEPEND' for %s scale distribution?"%sn_or_host)
    else:
        dist_dict['SCALE']=lambda z:scale_func()[0]
    return(dist_dict)