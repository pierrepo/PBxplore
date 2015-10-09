#!/usr/bin/env python

import json
import sys

data = json.load(sys.stdin)

fig_count = 0
for stuff in data[1]:
    if u't' in stuff and stuff[u't'] == 'Para':
        for substuff in  stuff[u'c']:
            if substuff[u't'] == 'Image':
                fig_count += 1
                if substuff[u'c'][0][0][u'c'] == u'png':
                    substuff[u'c'][0][0][u'c'] = u'Figure {}'.format(fig_count)
                if substuff[u'c'][1][0].startswith(u'./notebooks/'):
                    substuff[u'c'][1][0] = './' + substuff[u'c'][1][0][12:]

json.dump(data, sys.stdout)
