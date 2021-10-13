import pyVPM
import unittest

class TestPoint2d(unittest.TestCase):
    def test_simple_point(self):
        point = pyVPM.Point2d(0.5, 0.4)
        self.assertEqual(0.5, point.x)
        self.assertEqual(0.4, point.y)

if __name__ == '__main__':
    unittest.main()
