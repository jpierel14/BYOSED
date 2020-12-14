import sys,os
_SRCDIR_ = os.path.abspath(os.path.join(
	os.path.dirname(os.path.abspath(__file__)),'..'))
sys.path.append(_SRCDIR_)
import byosed

def test_byosed():
	mySED=byosed.GeneralSED()
	flux = mySED.warp_SED()
	mod = mySED.to_sn_model()

if __name__ == '__main__':
	test_byosed()