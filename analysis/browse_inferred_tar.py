import cv2 as cv
import sys
import csv
import numpy as np
import math
import tarfile
import os
import glob
import time

# 224 hardcoded in ryby-infer scripts
oc=224.0
# show original image enlarged by
scaling_factor = 2.0

img_dims = ( round(oc*scaling_factor), round(oc*scaling_factor) )

# graph on bottom of the window
graph_height = 56.0

canvas = np.zeros( (round((oc+graph_height)*scaling_factor), round(oc*scaling_factor),3), np.uint8)

looping = True
RL = False
low = -10
high = 10

def del_browse_running():
    if os.path.exists("browse_running"):
        print("deleting browse_running")
        os.remove("browse_running")

try:
    records = []
    displacements = []
    with open(glob.glob(sys.argv[1])[0]) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        skip = True
        for row in reader:
          if skip:
            skip = False
            continue
          records.append(row)

    tar0 = tarfile.open(sys.argv[2])
    cmdline_dir = sys.argv[2].split("/")[0]
    cmdline_fishnr = sys.argv[2].split("/")[1]

    with open(sys.argv[3]) as framefile:
        try:
            frame_arr = framefile.read().split(";")
            gotoframe = int(frame_arr[0])
            if (len(frame_arr) == 5):
                framenr_dir = frame_arr[1]
                framenr_fishnr = frame_arr[2].strip()
                low = int(frame_arr[3])
                high = int(frame_arr[4])
                print("Showing frame", gotoframe, "for fish#", framenr_fishnr, "from dataset", framenr_dir)
            else:
                tar0.close()
                del_browse_running()
                sys.exit("Quitting, error parsing framenr.csv file.")
        except SystemExit:
            raise
        except:
            print("Fallback to frame 1")
            gotoframe = 1

    slider_max = len(records)-1
    binsize = int(slider_max/round(oc*scaling_factor))

    binnr=0
    for idx,row in enumerate(records):
        if (idx == 0):
            dispmax = 0.0
            continue
        if (idx // binsize == binnr):
            disp=math.sqrt( (float(records[idx][2])-float(records[idx-1][2]))**2.0 + (float(records[idx][3])-float(records[idx-1][3]))**2.0 )
            if (disp > dispmax):
                dispmax = disp
        else:
            displacements.append(int(dispmax))
            binnr +=1
            dispmax = 0.0
    oamin = min(displacements)
    oamax = max(displacements)
    displacements = [ int(graph_height*scaling_factor*(oamax-x)/(oamax-oamin) ) for x in displacements ]
except:
    print("\nRun as: " + sys.argv[0] + " some_inferred.csv images.tar framenr.csv\n")
    raise

window_title = 'Inferred head/tail positions'

def get_np_array_from_tar_object(tar_extractfl):
    return np.asarray(bytearray(tar_extractfl.read()), dtype=np.uint8)

def on_trackbar(frame):
    try:
        img = None
        img = cv.imdecode(get_np_array_from_tar_object(tar0.extractfile(records[frame][1])), 3)
        if img is None:
            sys.exit("Could not read the image.")
        scaled_img = cv.resize(img, img_dims)
        scaled_img.resize(round((oc+graph_height)*scaling_factor),round(oc*scaling_factor),3)
        cv.drawMarker(scaled_img,(round(float(records[frame][2])*scaling_factor),round(float(records[frame][3])*scaling_factor)), (0,255,255), cv.MARKER_CROSS, round(8*scaling_factor), 2)
        cv.drawMarker(scaled_img,(round(float(records[frame][4])*scaling_factor),round(float(records[frame][5])*scaling_factor)), (255,255,0), cv.MARKER_CROSS, round(8*scaling_factor), 2)
        for idx,disp in enumerate(displacements):
            if (idx == 0):
                continue
            if (idx+1 < round(oc*scaling_factor)):
                cv.line(scaled_img,(idx, round(oc*scaling_factor)+disp),(idx-1, round(oc*scaling_factor)+displacements[idx-1]),(255,255,255),1)
#                scaled_img[round(224.0*scaling_factor)+disp,idx]=(0,0,255)
        binnr = frame // binsize
        if (binnr < round(oc*scaling_factor)):
            cv.circle(scaled_img,( binnr , round(oc*scaling_factor)+displacements[binnr] ), 3, (0,0,255), -1)
        cv.imshow(window_title, scaled_img)
    except OSError:
        sys.exit("Could not read the image.")
    except SystemExit:
        raise
    except:
        raise

def run_loop(*args):
    global looping
    looping = not looping
    print("Looping set to", looping)
    if (not looping):
        gotoframe = 1

def switch_RL(*args):
    global RL
    RL = not RL

def set_around(*args):
    global low
    global high
    if (args[1] == "decrement"):
        if (RL):
            high -= 5
            if (high < 0):
                high = 0
        else:
            low += 5
            if (low > 0):
                low = 0
    if (args[1] == "increment"):
        if (RL):
            high += 5
        else:
            low -= 5
    print("Showing", -1*low, "frames before and" , high, "frames after")

try:
    cv.namedWindow(window_title)

    trackbar_name = 'Frames'
    cv.createTrackbar(trackbar_name, window_title , 1, slider_max+1, on_trackbar)

    #print("Jump to frame",gotoframe)
    on_trackbar(gotoframe-1)

    cv.createButton("loop",run_loop,None,cv.QT_CHECKBOX|cv.QT_NEW_BUTTONBAR,True)
    cv.createButton("<<",set_around,"decrement",cv.QT_PUSH_BUTTON|cv.QT_NEW_BUTTONBAR,False)
    cv.createButton("R/L",switch_RL,None,cv.QT_CHECKBOX,True)
    cv.createButton(">>",set_around,"increment",cv.QT_PUSH_BUTTON,False)
except:
    del_browse_running()

try:
    filename = "browse_running"
    if not os.path.exists(filename):
        open(filename, 'w').close()
except:
    tar0.close()
    sys.exit("Could not create \"browse_running\" file.")


while cv.getWindowProperty(window_title,cv.WND_PROP_VISIBLE) > 0:
    if (gotoframe != 1):
        cv.setTrackbarPos(trackbar_name, window_title, gotoframe-1)
    if (looping):
        with open(sys.argv[3]) as framefile:
            try:
                frame_arr = framefile.read().split(";")
                newframe = int(frame_arr[0])
                if (len(frame_arr) == 5):
                    framenr_dir = frame_arr[1]
                    framenr_fishnr = frame_arr[2].strip()
                    low = int(frame_arr[3])
                    if (newframe+low < 0):
                      low = 0
                    high = int(frame_arr[4])
                    if (framenr_dir != cmdline_dir):
                        tar0.close()
                        del_browse_running()
                        sys.exit("Quitting, directory changed,", framenr_dir, "!=", cmdline_dir)
                    if (framenr_fishnr != cmdline_fishnr):
                        tar0.close()
                        del_browse_running()
                        sys.exit("Quitting, fish number changed,", framenr_fishnr, "!=", cmdline_fishnr)
            except SystemExit:
                raise
            except:
                newframe = 1
    if (gotoframe != 1):
        if (looping):
            for fr in range(newframe+low, newframe+high):
                cv.setTrackbarPos(trackbar_name, window_title, fr-1)
                k = cv.waitKey(100)
                if k == ord("q"):
                    tar0.close()
                    del_browse_running()
                    sys.exit("Quitting, \"q\" pressed.")
                if (not looping):
                    break
            gotoframe = newframe

    k = cv.waitKey(200)

    if k == ord("q"):
        tar0.close()
        del_browse_running()
        sys.exit("Quitting, \"q\" pressed.")
else:
    tar0.close()
    del_browse_running()
    sys.exit("Window closed, quitting.")
