#include "Structure_Wing.hpp"
#include <iostream>

namespace VPM
{

    Structure_Wing::Structure_Wing()
    {
       m_origo = Point2d(.45, .0);
       m_charlength = .12;
       m_p.clear();
       m_p.push_back(Point2d( 1.00000000,  0.00000000+0.0128347)); 
       m_p.push_back(Point2d( 0.99939800,  0.00003500+0.0128347));
       m_p.push_back(Point2d( 0.99759200,  0.00013700+0.0128347));
       m_p.push_back(Point2d( 0.99458800,  0.00029600+0.0128347));
       m_p.push_back(Point2d( 0.99039300,  0.00049700+0.0128347));
       m_p.push_back(Point2d( 0.98501600,  0.00071900+0.0128347));
       m_p.push_back(Point2d( 0.97847000,  0.00093500+0.0128347));
       m_p.push_back(Point2d( 0.97077200,  0.00111200+0.0128347));
       m_p.push_back(Point2d( 0.96194000,  0.00121200+0.0128347));
       m_p.push_back(Point2d( 0.95199500,  0.00119700+0.0128347));
       m_p.push_back(Point2d( 0.94096100,  0.00103300+0.0128347));
       m_p.push_back(Point2d( 0.92886400,  0.00069400+0.0128347));
       m_p.push_back(Point2d( 0.91573500,  0.00015700+0.0128347));
       m_p.push_back(Point2d( 0.90160400, -0.00060000+0.0128347));
       m_p.push_back(Point2d( 0.88650500, -0.00159200+0.0128347));
       m_p.push_back(Point2d( 0.87047600, -0.00282900+0.0128347));
       m_p.push_back(Point2d( 0.85355300, -0.00431400+0.0128347));
       m_p.push_back(Point2d( 0.83577900, -0.00604800+0.0128347));
       m_p.push_back(Point2d( 0.81719700, -0.00802700+0.0128347));
       m_p.push_back(Point2d( 0.79785000, -0.01024400+0.0128347));
       m_p.push_back(Point2d( 0.77778500, -0.01269000+0.0128347));
       m_p.push_back(Point2d( 0.75705100, -0.01535700+0.0128347));
       m_p.push_back(Point2d( 0.73569800, -0.01823200+0.0128347));
       m_p.push_back(Point2d( 0.71377800, -0.02128900+0.0128347));
       m_p.push_back(Point2d( 0.69134200, -0.02449500+0.0128347));
       m_p.push_back(Point2d( 0.66844500, -0.02781400+0.0128347));
       m_p.push_back(Point2d( 0.64514200, -0.03120700+0.0128347));
       m_p.push_back(Point2d( 0.62149000, -0.03463100+0.0128347));
       m_p.push_back(Point2d( 0.59754500, -0.03804300+0.0128347));
       m_p.push_back(Point2d( 0.57336500, -0.04139700+0.0128347));
       m_p.push_back(Point2d( 0.54900900, -0.04464200+0.0128347));
       m_p.push_back(Point2d( 0.52453400, -0.04771900+0.0128347));
       m_p.push_back(Point2d( 0.50000000, -0.05056300+0.0128347));
       m_p.push_back(Point2d( 0.47546600, -0.05309900+0.0128347));
       m_p.push_back(Point2d( 0.45099100, -0.05525700+0.0128347));
       m_p.push_back(Point2d( 0.42663500, -0.05697900+0.0128347));
       m_p.push_back(Point2d( 0.40245500, -0.05822400+0.0128347));
       m_p.push_back(Point2d( 0.37851000, -0.05897400+0.0128347));
       m_p.push_back(Point2d( 0.35485800, -0.05923600+0.0128347));
       m_p.push_back(Point2d( 0.33155500, -0.05904600+0.0128347));
       m_p.push_back(Point2d( 0.30865800, -0.05845900+0.0128347));
       m_p.push_back(Point2d( 0.28622200, -0.05754700+0.0128347));
       m_p.push_back(Point2d( 0.26430200, -0.05637600+0.0128347));
       m_p.push_back(Point2d( 0.24294900, -0.05499400+0.0128347));
       m_p.push_back(Point2d( 0.22221500, -0.05342700+0.0128347));
       m_p.push_back(Point2d( 0.20215000, -0.05169400+0.0128347));
       m_p.push_back(Point2d( 0.18280300, -0.04980500+0.0128347));
       m_p.push_back(Point2d( 0.16422100, -0.04777300+0.0128347));
       m_p.push_back(Point2d( 0.14644700, -0.04561000+0.0128347));
       m_p.push_back(Point2d( 0.12952400, -0.04332600+0.0128347));
       m_p.push_back(Point2d( 0.11349500, -0.04092900+0.0128347));
       m_p.push_back(Point2d( 0.09839600, -0.03843100+0.0128347));
       m_p.push_back(Point2d( 0.08426500, -0.03584300+0.0128347));
       m_p.push_back(Point2d( 0.07113600, -0.03317000+0.0128347));
       m_p.push_back(Point2d( 0.05903900, -0.03041600+0.0128347));
       m_p.push_back(Point2d( 0.04800500, -0.02758600+0.0128347));
       m_p.push_back(Point2d( 0.03806000, -0.02468500+0.0128347));
       m_p.push_back(Point2d( 0.02922800, -0.02172200+0.0128347));
       m_p.push_back(Point2d( 0.02153000, -0.01870700+0.0128347));
       m_p.push_back(Point2d( 0.01498400, -0.01564900+0.0128347));
       m_p.push_back(Point2d( 0.00960700, -0.01255900+0.0128347));
       m_p.push_back(Point2d( 0.00541200, -0.00944300+0.0128347));
       m_p.push_back(Point2d( 0.00240800, -0.00630800+0.0128347));
       m_p.push_back(Point2d( 0.00060200, -0.00316000+0.0128347));
       m_p.push_back(Point2d( 0.00000000,  0.00000000+0.0128347));
       m_p.push_back(Point2d( 0.00060200,  0.00316500+0.0128347));
       m_p.push_back(Point2d( 0.00240800,  0.00630600+0.0128347));
       m_p.push_back(Point2d( 0.00541200,  0.00941600+0.0128347));
       m_p.push_back(Point2d( 0.00960700,  0.01248000+0.0128347));
       m_p.push_back(Point2d( 0.01498400,  0.01548900+0.0128347));
       m_p.push_back(Point2d( 0.02153000,  0.01844100+0.0128347));
       m_p.push_back(Point2d( 0.02922800,  0.02134800+0.0128347));
       m_p.push_back(Point2d( 0.03806000,  0.02421900+0.0128347));
       m_p.push_back(Point2d( 0.04800500,  0.02706200+0.0128347));
       m_p.push_back(Point2d( 0.05903900,  0.02987400+0.0128347));
       m_p.push_back(Point2d( 0.07113600,  0.03264400+0.0128347));
       m_p.push_back(Point2d( 0.08426500,  0.03536000+0.0128347));
       m_p.push_back(Point2d( 0.09839600,  0.03801100+0.0128347));
       m_p.push_back(Point2d( 0.11349500,  0.04058500+0.0128347));
       m_p.push_back(Point2d( 0.12952400,  0.04307100+0.0128347));
       m_p.push_back(Point2d( 0.14644700,  0.04545700+0.0128347));
       m_p.push_back(Point2d( 0.16422100,  0.04772900+0.0128347));
       m_p.push_back(Point2d( 0.18280300,  0.04987400+0.0128347));
       m_p.push_back(Point2d( 0.20215000,  0.05188500+0.0128347));
       m_p.push_back(Point2d( 0.22221500,  0.05375300+0.0128347));
       m_p.push_back(Point2d( 0.24294900,  0.05547000+0.0128347));
       m_p.push_back(Point2d( 0.26430200,  0.05702600+0.0128347));
       m_p.push_back(Point2d( 0.28622200,  0.05841400+0.0128347));
       m_p.push_back(Point2d( 0.30865800,  0.05962900+0.0128347));
       m_p.push_back(Point2d( 0.33155500,  0.06066000+0.0128347));
       m_p.push_back(Point2d( 0.35485800,  0.06149700+0.0128347));
       m_p.push_back(Point2d( 0.37851000,  0.06213300+0.0128347));
       m_p.push_back(Point2d( 0.40245500,  0.06256200+0.0128347));
       m_p.push_back(Point2d( 0.42663500,  0.06277900+0.0128347));
       m_p.push_back(Point2d( 0.45099100,  0.06277400+0.0128347));
       m_p.push_back(Point2d( 0.47546600,  0.06253000+0.0128347));
       m_p.push_back(Point2d( 0.50000000,  0.06202900+0.0128347));
       m_p.push_back(Point2d( 0.52453400,  0.06125400+0.0128347));
       m_p.push_back(Point2d( 0.54900900,  0.06019400+0.0128347));
       m_p.push_back(Point2d( 0.57336500,  0.05884500+0.0128347));
       m_p.push_back(Point2d( 0.59754500,  0.05721800+0.0128347));
       m_p.push_back(Point2d( 0.62149000,  0.05534400+0.0128347));
       m_p.push_back(Point2d( 0.64514200,  0.05325800+0.0128347));
       m_p.push_back(Point2d( 0.66844500,  0.05099300+0.0128347));
       m_p.push_back(Point2d( 0.69134200,  0.04857500+0.0128347));
       m_p.push_back(Point2d( 0.71377800,  0.04602900+0.0128347));
       m_p.push_back(Point2d( 0.73569800,  0.04337700+0.0128347));
       m_p.push_back(Point2d( 0.75705100,  0.04064100+0.0128347));
       m_p.push_back(Point2d( 0.77778500,  0.03784700+0.0128347));
       m_p.push_back(Point2d( 0.79785000,  0.03501700+0.0128347));
       m_p.push_back(Point2d( 0.81719700,  0.03217600+0.0128347));
       m_p.push_back(Point2d( 0.83577900,  0.02934700+0.0128347));
       m_p.push_back(Point2d( 0.85355300,  0.02655400+0.0128347));
       m_p.push_back(Point2d( 0.87047600,  0.02381700+0.0128347));
       m_p.push_back(Point2d( 0.88650500,  0.02115300+0.0128347));
       m_p.push_back(Point2d( 0.90160400,  0.01858000+0.0128347));
       m_p.push_back(Point2d( 0.91573500,  0.01611300+0.0128347));
       m_p.push_back(Point2d( 0.92886400,  0.01376900+0.0128347));
       m_p.push_back(Point2d( 0.94096100,  0.01156200+0.0128347));
       m_p.push_back(Point2d( 0.95199500,  0.00950800+0.0128347));
       m_p.push_back(Point2d( 0.96194000,  0.00762200+0.0128347));
       m_p.push_back(Point2d( 0.97077200,  0.00591500+0.0128347));
       m_p.push_back(Point2d( 0.97847000,  0.00440100+0.0128347));
       m_p.push_back(Point2d( 0.98501600,  0.00309200+0.0128347));
       m_p.push_back(Point2d( 0.99039300,  0.00200100+0.0128347));
       m_p.push_back(Point2d( 0.99458800,  0.00113700+0.0128347));
       m_p.push_back(Point2d( 0.99759200,  0.00051000+0.0128347));
       m_p.push_back(Point2d( 0.99939800,  0.00012800+0.0128347));
       m_p.push_back(Point2d( 1.00000000,  0.00000000+0.0128347));

    }
    Structure_Wing::~Structure_Wing()
    {
    }

    namespace {
        //struct Segment {
        //    Point2d p1;
        //    Point2d p2;
        //};

        //bool onSegment(Segment l1, Point2d p)
        //{   //check whether p is on the Segment or not
        //    bool ret = false;
        //    if ( p.x <= std::max(l1.p1.x, l1.p2.x) &&
        //         p.x <= std::min(l1.p1.x, l1.p2.x) &&
        //        (p.y <= std::max(l1.p1.y, l1.p2.y) && p.y <= std::min(l1.p1.y, l1.p2.y))
        //        )
        //    {
        //       ret = true;
        //    }
        //    return ret;
        //}

        //int direction(Point2d a, Point2d b, Point2d c)
        //{
        //    int ret;
        //    int val = (b.y-a.y)*(c.x-b.x)-(b.x-a.x)*(c.y-b.y);
        //    if (val == 0)
        //    {
        //        ret = 0;     //colinear
        //    }
        //    else if(val < 0)
        //    {
        //        ret = 2;    //anti-clockwise direction
        //    }
        //    else
        //    {
        //        ret = 1;    //clockwise direction
        //    }
        //    return ret;
        //}

        //bool isIntersect(Segment l1, Segment l2)
        //{
        //    bool ret = false;
        //    //four direction for two Segments and points of other Segment
        //    int dir1 = direction(l1.p1, l1.p2, l2.p1);
        //    int dir2 = direction(l1.p1, l1.p2, l2.p2);
        //    int dir3 = direction(l2.p1, l2.p2, l1.p1);
        //    int dir4 = direction(l2.p1, l2.p2, l1.p2);

        //    if(dir1 != dir2 && dir3 != dir4)
        //    {
        //        ret = true; //they are intersecting
        //    }
        //    else if(dir1==0 && onSegment(l1, l2.p1))
        //    { //when p2 of Segment2 is on Segment1
        //        ret = true;
        //    }
        //    else if(dir2==0 && onSegment(l1, l2.p2)) //when p1 of Segment2 is on Segment1
        //    {
        //        ret = true;
        //    }
        //    else if(dir3==0 && onSegment(l2, l1.p1)) //when p2 of Segment1 is on Segment2
        //    {
        //        ret = true;
        //    }
        //    else if(dir4==0 && onSegment(l2, l1.p2)) //when p1 of Segment1 is on Segment2
        //    {
        //        ret = true;
        //    }
        //    return ret;
        //}
        bool ccw(Point2d A, Point2d B, Point2d C)
        {
            return (C.y-A.y) * (B.x-A.x) > (B.y-A.y) * (C.x-A.x);
        }
        bool intersect(Point2d A, Point2d B, Point2d C, Point2d D)
        {
            return (ccw(A,C,D) != ccw(B,C,D)) && (ccw(A,B,C) != ccw(A,B,D));
        }
    }

    bool Structure_Wing::isInside(const Point2d pos, const double pad)
    {

       bool ret = false;
           int count_intersections = 0;
       if (L2Norm(pos-Point2d(0,0))<1e-12)
       {
           ret = true;
       }
       else
       {
           Point2d p1 = pos;
           Point2d p2 = Point2d(pos.x, +0.1);
           //Point2d p1 = Point2d(-0.1, pos.y);
           //Point2d p2 = Point2d(+1.1, pos.y);
           for ( int i = 0; i < m_p.size()-1; i++)
           {
               Point2d s1 = m_p[i];
               Point2d s2 = m_p[i+1];
               if (intersect(p1, p2, s1, s2))
               //if (isIntersect(s1, s2))
               {
                   count_intersections +=1;
               }
           }
           if (count_intersections%2 != 0)
           {
               ret = true;
           }
       }
       //std::cerr<<"P=("<<pos.x<<","<<pos.y<<")="<<count_intersections<<std::endl;
       return ret;
    }

    void Structure_Wing::getOrigo(Point2d & pos)
    {
        pos = m_origo;
    }

    double Structure_Wing::getCharacteristicLength()
    {
        return m_charlength;
    }


}
