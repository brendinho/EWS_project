#!/bin/bash

sshpass -p "Brendon1993" ssh -tt -o StrictHostKeyChecking=no b2philli@graham.sharcnet.ca << EOF
        cd /home/b2philli;
	# touch logged_in.txt;
	. release_the_kraken_SHARCNET.sh;
	exit;
EOF

