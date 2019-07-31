##Tests for pipeline
import sys,os
sys.path.append(os.path.join(os.path.abspath(__file__),'..','..'))
import byosed

def test_byosed():
	mySED=byosed.GeneralSED()
	flux = mySED.warp_SED()
	mod = mySED.to_sn_model()

