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

    dist = skewed_normal(name,a=a,b=b)
    sample=np.arange(a,b,.01)
    return(lambda : np.random.choice(sample,1,
                                     p=dist._pdf(sample,dist_dat[dist_type+'_DIST_PEAK'],dist_dat[dist_type+'_DIST_SIGMA'][0],dist_dat[dist_type+'_DIST_SIGMA'][1])))


def _param_from_dist(dist_file,path):
    dist=np.loadtxt(os.path.join(path,dist_file))
    a=np.min(dist)-abs(np.min(dist))
    b=np.max(dist)+abs(np.max(dist))
    sample=np.arange(a,b,.01)
    pdf=gaussian_kde(dist.T).pdf(np.arange(a,b,.01))
    return(lambda : np.random.choice(sample,1,p=pdf/np.sum(pdf)))

def get_distribution(name,dist_dat,path,sn_or_host):
    dist_dict={}
    if np.any(['DIST_FILE' in x for x in dist_dat.keys() if 'PARAM' in x]):
        if sn_or_host+'_PARAM_DIST_FILE' in dist_dat.keys():
            dist_dict['PARAM']=_param_from_dist(dist_dat[sn_or_host+'_PARAM_DIST_FILE'],path)
        else:
            raise RuntimeError("You may have a typo, did you mean to set 'PARAM_DIST_FILE' for %s?"%sn_or_host)


    elif np.any(['PARAM' in x for x in dist_dat.keys()]):
        try:
            dist_dict['PARAM']=_skewed_normal(name,dist_dat,sn_or_host+'_PARAM')
        except:
            raise RuntimeError("You may have a typo in the variables of your %s param distribution(s)."%sn_or_host)

    if np.any(['DIST_FILE' in x for x in dist_dat.keys() if 'SCALE' in x]):
        if 'SCALE_DIST_FILE' in dist_dat.keys():
            dist_dict['SCALE']=_param_from_dist(dist_dat['SCALE_DIST_FILE'],path)
        else:
            raise RuntimeError("You may have a typo, did you mean to set 'SCALE_DIST_FILE'?")

    elif np.any(['SCALE' in x for x in dist_dat.keys()]):
        try:
            dist_dict['SCALE']=_skewed_normal(name,dist_dat,'SCALE')
        except:
            raise RuntimeError("You may have a typo in the variables of your %s scale distribution(s)."%sn_or_host)
    else:
        raise RuntimeError("Must supply scale distribution for every effect.")

    return(dist_dict)