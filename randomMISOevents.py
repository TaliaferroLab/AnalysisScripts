#Takes MISO bayes file output and picks random events from that file.  Returns list of random events

#Usage: python randomMISOevents <file.miso_bf> <number_of_random_events> <output.txt>

import sys
import random

def randomMISOevents(eventfile, number):
    number = int(number)
    events = []
    eventsfh = open(eventfile, 'r')
    for event in eventsfh:
        event = event.strip().split('\t')
        if event[0] != 'event_name':
            events.append(event[0])
        else:
            print 'Skipping header line'

    eventsfh.close()

    print 'Have list of %i events' % len(events)
    print 'Choosing %i events at random' % number

    randomevents = random.sample(events, number)

    return randomevents

outfh = open(sys.argv[3], 'w')
for event in randomMISOevents(sys.argv[1], sys.argv[2]):
    outfh.write(str(event) + '\n')

outfh.close()
