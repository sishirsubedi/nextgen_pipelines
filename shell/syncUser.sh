rsync -avrP /etc/passwd root@n001:/etc
rsync -avrP /etc/group root@n001:/etc
rsync -avrP /etc/shadow root@n001:/etc

rsync -avrP /etc/passwd root@n002:/etc
rsync -avrP /etc/group root@n002:/etc
rsync -avrP /etc/shadow root@n002:/etc

rsync -avrP /etc/passwd root@n003:/etc
rsync -avrP /etc/group root@n003:/etc
rsync -avrP /etc/shadow root@n003:/etc

rsync -avrP /etc/passwd root@n004:/etc
rsync -avrP /etc/group root@n004:/etc
rsync -avrP /etc/shadow root@n004:/etc

ssh root@n001 nfsidmap -c
ssh root@n001 /etc/init.d/rpcidmapd restart

ssh root@n002 nfsidmap -c
ssh root@n002 /etc/init.d/rpcidmapd restart

ssh root@n002 nfsidmap -c
ssh root@n002 /etc/init.d/rpcidmapd restart

ssh root@n002 nfsidmap -c
ssh root@n002 /etc/init.d/rpcidmapd restart

