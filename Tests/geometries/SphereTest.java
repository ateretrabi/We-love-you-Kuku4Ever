package geometries;

import org.junit.jupiter.api.Test;
import primitives.Point3D;
import primitives.Ray;
import primitives.Vector;

import java.util.LinkedList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * @author Dan zilberstein
 */
class SphereTest {

    Sphere sphere = new Sphere(1d, new Point3D(1, 0, 0));


    @Test
    void getNormalTest1() {
        Sphere sp = new Sphere(1.0, new Point3D(0, 0, 1));
        assertEquals(new Vector(0,0,1),sp.getNormal(new Point3D(0,0,2)));
    }

    @Test
    void getNormalTest2() {
        Sphere sp = new Sphere(1,new Point3D(0,0,1));
        assertNotEquals(new Vector(0,0,1),sp.getNormal(new Point3D(0,1,1)));
        System.out.println(sp.getNormal(new Point3D(0,1,1)));
    }

    @Test
    void findIntersectionsTestEP_1()
   {
        // ============ Equivalence Partitions Tests ==============
       Point3D p1 = new Point3D(0.0651530771650466, 0.355051025721682, 0);
       Point3D p2 = new Point3D(1.53484692283495, 0.844948974278318, 0);
       List<Point3D> exp = List.of(p1, p2);

       // TC01: Ray's line is outside the sphere (0 points)
       assertNull(sphere.findIntersections(new Ray(new Point3D(-1, 0, 0), new Vector(1, 1, 0))),
               "Ray's line out of sphere");

   }

    /**
     * @author Dan Zilberstein
     */
    @Test
    public void findIntersectionsTest() {
//        Sphere sphere = new Sphere(1d, new Point3D(1, 0, 0));

        // ============ Equivalence Partitions Tests ==============
        Point3D p1 = new Point3D(0.0651530771650466, 0.355051025721682, 0);
        Point3D p2 = new Point3D(1.53484692283495, 0.844948974278318, 0);
        List<Point3D> exp = List.of(p1, p2);
        List<Point3D> points = new LinkedList<>();

        // TC01: Ray's line is outside the sphere (0 points)
        assertNull(sphere.findIntersections(new Ray(new Point3D(-1, 0, 0), new Vector(1, 1, 0))),
                "Ray's line out of sphere");

        // TC02: Ray starts before and crosses the sphere (2 points)
        List<Intersectable.GeoPoint> result02 = sphere.findIntersections(new Ray(new Point3D(-1, 0, 0), new Vector(3, 1, 0)));

        assertEquals(2, result02.size(), "Wrong number of points");

        points.clear();
        for (Intersectable.GeoPoint geo : result02) {
            points.add(geo.getPoint());
        }

        if (points.get(0).get_x().get() > points.get(1).get_x().get()) {
            points = List.of(points.get(1), points.get(0));
        }
        assertEquals(exp, points, "Ray crosses sphere");

        // TC03: Ray starts inside the sphere (1 point)
        List<Intersectable.GeoPoint> result03 = sphere.findIntersections(new Ray(new Point3D(0.5, 0.5, 0), new Vector(3, 1, 0)));
        points.clear();
        for (Intersectable.GeoPoint geo : result03) {
            points.add(geo.getPoint());
        }
        assertEquals(List.of(p2), points,
                "Ray from inside sphere");
        // TC04: Ray starts after the sphere (0 points)
        assertNull(sphere.findIntersections(new Ray(new Point3D(2, 1, 0), new Vector(3, 1, 0))),
                "Sphere behind Ray");

        // =============== Boundary Values Tests ==================
        // **** Group: Ray's line crosses the sphere (but not the center)

        // TC11: Ray starts at sphere and goes inside (1 points)
        List<Intersectable.GeoPoint> result11 = sphere.findIntersections(new Ray(new Point3D(1, -1, 0), new Vector(1, 1, 0)));
        points.clear();
        for (Intersectable.GeoPoint geo : result11) {
            points.add(geo.getPoint());
        }

        assertEquals(List.of(new Point3D(2, 0, 0)), points, "Ray from sphere inside");

        // TC12: Ray starts at sphere and goes outside (0 points)
        assertNull(sphere.findIntersections(new Ray(new Point3D(2, 0, 0), new Vector(1, 1, 0))),
                "Ray from sphere outside");

        // **** Group: Ray's line goes through the center
        // TC13: Ray starts before the sphere (2 points)
        List<Intersectable.GeoPoint> result13 = sphere.findIntersections(new Ray(new Point3D(1, -2, 0), new Vector(0, 1, 0)));

        assertEquals(2, result13.size(), "Wrong number of points");
        points.clear();
        for (Intersectable.GeoPoint geo : result13) {
            points.add(geo.getPoint());
        }
        if (points.get(0).get_y().get() > points.get(1).get_y().get()) {
            result13 = List.of(result13.get(1), result13.get(0));
        }
        assertEquals(
                List.of(new Point3D(1, -1, 0), new Point3D(1, 1, 0)),
                points,
                "Line through O, ray crosses sphere");

        // TC14: Ray starts at sphere and goes inside (1 points)
        List<Intersectable.GeoPoint> result14 = sphere.findIntersections(new Ray(new Point3D(1, -1, 0), new Vector(0, 1, 0)));
        points.clear();
        for (Intersectable.GeoPoint geo : result14) {
            points.add(geo.getPoint());
        }

        assertEquals(List.of(new Point3D(1, 1, 0)), points,
                "Line through O, ray from and crosses sphere");

        // TC15: Ray starts inside (1 points)
        List<Intersectable.GeoPoint> result15 = sphere.findIntersections(new Ray(new Point3D(1, 0.5, 0), new Vector(0, 1, 0)));
        points.clear();
        for (Intersectable.GeoPoint geo : result15) {
            points.add(geo.getPoint());
        }
        assertEquals(List.of(new Point3D(1, 1, 0)),
                points,
                "Line through O, ray from inside sphere");

        // TC16: Ray starts at the center (1 points)
        List<Intersectable.GeoPoint> result16 = sphere.findIntersections(new Ray(new Point3D(1, 0, 0), new Vector(0, 1, 0)));
        points.clear();
        for (Intersectable.GeoPoint geo : result16) {
            points.add(geo.getPoint());
        }
        assertEquals(List.of(new Point3D(1, 1, 0)), points,
                "Line through O, ray from O");

        // TC17: Ray starts at sphere and goes outside (0 points)
        assertNull(sphere.findIntersections(new Ray(new Point3D(1, 1, 0), new Vector(0, 1, 0))),
                "Line through O, ray from sphere outside");

        // TC18: Ray starts after sphere (0 points)
        assertNull(sphere.findIntersections(new Ray(new Point3D(1, 2, 0), new Vector(0, 1, 0))),
                "Line through O, ray outside sphere");

        // **** Group: Ray's line is tangent to the sphere (all tests 0 points)
        // TC19: Ray starts before the tangent point
        assertNull(sphere.findIntersections(new Ray(new Point3D(0, 1, 0), new Vector(1, 0, 0))),
                "Tangent line, ray before sphere");

        // TC20: Ray starts at the tangent point
        assertNull(sphere.findIntersections(new Ray(new Point3D(1, 1, 0), new Vector(1, 0, 0))),
                "Tangent line, ray at sphere");

        // TC21: Ray starts after the tangent point
        assertNull(sphere.findIntersections(new Ray(new Point3D(2, 1, 0), new Vector(1, 0, 0))),
                "Tangent line, ray after sphere");

        // **** Group: Special cases
        // TC19: Ray's line is outside, ray is orthogonal to ray start to sphere's
        // center line
        assertNull(sphere.findIntersections(new Ray(new Point3D(-1, 0, 0), new Vector(0, 0, 1))),
                "Ray orthogonal to ray head -> O line");

    }

}
