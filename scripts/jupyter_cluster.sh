## get tunneling info
XDG_RUNTIME_DIR=""
ipnport=$(shuf -i8000-9999 -n1)
ipnip=$(hostname -i)
port=$1

## print tunneling instructions to jupyter-log-{jobid}.txt
echo -e "
    Copy/Paste this in your local terminal to ssh tunnel with remote
    -----------------------------------------------------------------
    ssh -N -L $ipnport:$ipnip:$ipnport ivanov.vv@calc.cod.phystech.edu
    -----------------------------------------------------------------

    Then open a browser on your local machine to the following address
    ------------------------------------------------------------------
    localhost:$ipnport  (prefix w/ https:// if using password)
    ----------------------------------------------------------------- "-
    
## start an ipcluster instance and launch jupyter server
/usr/bin/ssh -N -f -R localhost:$ipnport:localhost:$ipnport ivanov.vv@calc.cod.phystech.edu
jupyter-notebook --no-browser --port=$ipnport --ip=$ipnip --notebook-dir='/home/common/ivanov.vv'

