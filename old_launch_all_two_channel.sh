#!/bin/bash

#sudo ./init_network.sh &

/home/radar/repos/usrp/stop_all.sh &

/bin/sleep 15

/home/radar/rst/bin/errlog -name kod.d -lp 41000 &
/home/radar/rst/bin/errlog -name kod.c -lp 43000 &

/home/radar/rst/bin/rawacfwrite -r kod.d -lp 41102 -ep 41000 &
/home/radar/rst/bin/rawacfwrite -r kod.c -lp 43102 -ep 43000 &

/home/radar/rst/bin/fitacfwrite -r kod.d -lp 41103 -ep 41000 &
/home/radar/rst/bin/fitacfwrite -r kod.c -lp 43103 -ep 43000 &


/home/radar/rst/bin/rtserver -rp 41104 -ep 41000 -tp 1401 & # ch 4
#/home/radar/rst/bin/rtserver -rp 43104 -ep 43000 -tp 1402 & # ch 3


# don't start remote watchdog on bootup
#if [ "$1" != "boot" ]
#then
ssh radar@kodiak-aux 'python3 /home/radar/SuperDARN_UHD_Server/tools/srr_watchdog.py &' &
#fi

/usr/bin/python3 /home/radar/repos/usrp/tools/srr_watchdog.py server &

/bin/sleep 30
/home/radar/rst/bin/schedule -l /data/ros/scdlog/kod.d.scdlog -f /data/ros/scd/kod.d.scd & 
/home/radar/rst/bin/schedule -l /data/ros/scdlog/kod.c.scdlog -f /data/ros/scd/kod.c.scd & 

#gdb -ex run --args uafscan --stid mcm -c 1 --debug # control program connects to usrp server
