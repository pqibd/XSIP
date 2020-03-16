sudo yum install -y https://centos7.iuscommunity.org/ius-release.rpm
sudo yum update
sudo yum -y install python36 python36-devel python36-pip python36-setuptools python36-tools python36-libs python36-tkinter
cd ~
python3 -m venv test_prj/test_env
source /home/pqi/test_prj/test_env/bin/activate
pip install -r test_prj/EDXAS/skes_requirements.txt

