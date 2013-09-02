#include "opencv2/opencv.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "circle_detection/tumult/math/Matrix.h"
#include "circle_detection/tumult/math/geom/ransac.h"
#include "ros/ros.h"
#include "std_msgs/String.h"
#include <circle_detection/circle_msg.h>

#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <csignal>

void signalHandler(int signum)
{
    exit(signum);
}


int lowThreshold = 25;
int houghThreshold = 50;

void DummyCallback(int, void*)
{
    if(lowThreshold == 0)
        lowThreshold = 1;
    if(houghThreshold == 0)
        houghThreshold = 1;
}

inline int max2(int a, int b)
{
    return ((a > b) ? a : b);
}

inline int min2(int a, int b)
{
    return ((a < b) ? a : b);
}

inline int max3(int a, int b, int c)
{
    return max2(max2(a,b),c);
}

inline int min3(int a, int b, int c)
{
    return min2(min2(a,b),c);
}

int main(int argc, char** argv)
{
    int ratio = 3;
    int kernel_size = 3;
    int sigma = 200;

    //connection to camera 
    CvCapture* capture = cvCaptureFromCAM(CV_CAP_ANY);

    if(!capture)
    {
        std::cerr << "ERROR: capture is NULL\n";
        getchar();
        return -1;
    }

    /* ONLY NEEDED FOR TESTING PURPOSES

       cvNamedWindow("Camera Image", CV_WINDOW_AUTOSIZE);
    // Create a Trackbar for user to enter threshold
    cv::createTrackbar("Min Threshold:","Camera Image" , &lowThreshold, max_lowThreshold, DummyCallback);
    // Create a Trackbar for user to enter threshold
    cv::createTrackbar("Hough Threshold:","Camera Image" , &houghThreshold, max_houghThreshold, DummyCallback);
    // Create a Trackbar for user to enter threshold
    cv::createTrackbar("Sigma:","Camera Image" , &sigma, 1000, DummyCallback);

    cvNamedWindow("Edge Image", CV_WINDOW_AUTOSIZE);

*/

    ros::init(argc,argv, "circle_publisher");
    ros::NodeHandle nh;

    //Let's publish messages as defined in /msg/cirlce.msg on a topic called 'detected_circles' with a max. buffer of 1000 messages 
    ros::Publisher circle_publisher = nh.advertise<circle_detection::circle_msg>("detected_circles",1000); 

    ros::Rate loop_rate(10);

    //used to create an Image Id
    size_t id_count = 0;

    while(ros::ok)
    {
        // Get a frame
        IplImage* frame = cvQueryFrame(capture);

        if (!frame)
        {
            std::cerr << "ERROR: frame is null...\n" ;
            getchar();
            break;
        }

        // image id 
        id_count++;

        cv::Mat src(frame);
        cv::Mat src_gray;

        cv::Mat dst;
        cv::Mat detected_edges;

        dst.create(src.size(), src.type());

        // covert the image to gray
        cvtColor(src,src_gray,CV_BGR2GRAY);

        // Reduce the noise so we avoid false circle detection
        GaussianBlur(src_gray, detected_edges, cv::Size(9, 9), sigma/100.0, 0);

        equalizeHist(detected_edges, detected_edges);

        // Canny detector
        Canny(detected_edges, detected_edges, lowThreshold, lowThreshold*ratio, kernel_size);

        // Using Canny's output as a mask, we display our result
        dst = cv::Scalar::all(0);

        src.copyTo(dst, detected_edges);

        std::vector<Point> edgePoints;


        // iterate through the  pixels of the canny image
        for(int j=0;j<detected_edges.cols;j++)
        {
            for(int i=0;i<detected_edges.rows;i++)
            {
                unsigned char &color = *(detected_edges.data+detected_edges.step*i + j*detected_edges.elemSize());
                unsigned char &bb = *(src.data+src.step*i + j*src.elemSize());
                unsigned char &gg = *(src.data+src.step*i + j*src.elemSize() + 1);
                unsigned char &rr = *(src.data+src.step*i + j*src.elemSize() + 2);

                // check if the pixel is black or white (only edges are white)
                if(color)
                {
                    int max = max3(rr,gg,bb);
                    int min = min3(rr,gg,bb);
                    int delta = max - min;

                    // check saturation (only colorfull circles will be detacted
                    if(delta > 20)
                    {
                        edgePoints.push_back(Point((double) j,(double) i));
                    }
                    else
                    {
                        // mark pixel as no longer relevant, i.e. the pixel isn't recognized as edge
                        color = 0;
                    }
                }
            }
        }


        std::vector<Circle> detectedCircles;

        /* Apply the RANSAC algorithm to find circles in the camera image
         * Paramters: sink vector, points, eps, iterations, minSupporters, maxCircles
         */ 
        ransac(detectedCircles, edgePoints, 5.0, 10000, 100, 1);

        /* ONLY NEDDED FOR TESTIGN PURPOSES

        // Draw the circles detected
        for(size_t i = 0; i < detectedCircles.size(); i++)
        {
        cv::Point center(cvRound(detectedCircles[i].center.x), cvRound(detectedCircles[i].center.y));
        int radius = cvRound(detectedCircles[i].radius);

        //            assert(radius > 0);

        // circle center
        circle(src, center, 3, cv::Scalar(0,255,0), -1, 8, 0 );
        // circle outline
        circle(src, center, radius, cv::Scalar(0,0,255), 3, 8, 0 );
        }

        //Results
        imshow("Camera Image", src);
        imshow("Edge Image", detected_edges);

*/


        for(size_t i=0;i < detectedCircles.size();i++)
        {
            circle_detection::circle_msg msg;
            std::stringstream ss;

            ss << "IMG" << id_count;

            msg.image_id = ss.str();

            unsigned int radius = cvRound(detectedCircles[i].radius);

            msg.radius = radius;

            msg.center_x = cvRound(detectedCircles[i].center.x) ;
            msg.center_y = cvRound(detectedCircles[i].center.y);
            circle_publisher.publish(msg);

            ros::spinOnce();
            loop_rate.sleep();

        }



        //Do not release the frame!
        //If ESC key pressed, Key=0x10001B under OpenCV 0.9.7(linux version),
        //remove higher bits using AND operator
        //cvWaitKey() is used as delay between the frames
        if ( (cvWaitKey(100) & 255) == 27 ) break;
    }

    /* ONLY NEEDED FOR TESTING PURPOSES
     
    // Release the capture device housekeeping
       cvReleaseCapture( &capture );
       cvDestroyWindow( "Display Image" );

    */

    return 0;
}
