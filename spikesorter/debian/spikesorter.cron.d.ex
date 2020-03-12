#
# Regular cron jobs for the spikesorter package
#
0 4	* * *	root	[ -x /usr/bin/spikesorter_maintenance ] && /usr/bin/spikesorter_maintenance
