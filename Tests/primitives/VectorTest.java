package primitives;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;
import static primitives.Util.isZero;

class VectorTest {

    @Test
    void testadd() {
        Vector v1 = new Vector(1.0, 1.0, 1.0);
        Vector v2 = new Vector(-1.0, -1.0, -1.5);

        v1 = v1.add(v2);
        assertEquals(new Vector(0.0,0.0,-0.5),v1);

        v2 = v2.add(v1);
        assertEquals(new Vector(-1.0, -1.0, -2.0),v2);

    }

    @Test
    void testsubtract() {
    }

    @Test
    void testscale() {
        Vector v1 = new Vector(1,1,1);

        Vector v1test = v1.scale(1);
        assertEquals(new Vector(1,1,1),v1test);

        Vector v2test = v1test.scale(2);
        assertEquals(new Vector(2,2,2),v2test);

        Vector v3test = v2test.scale(-2);
        assertEquals(new Vector(-4,-4,-4),v3test);
    }

    @Test
    void dotProduct() {
    }

    @Test
    void testcrossProduct() {

        Vector v1 = new Vector(3.5, -5.0, 10.0);
        Vector v2 = new Vector(2.5,7,0.5);
        Vector v3 = v1.crossProduct(v2);

        assertEquals( 0, v3.dotProduct(v2), 1e-10);
        assertEquals( 0, v3.dotProduct(v1), 1e-10);

        Vector v4 = v2.crossProduct(v1);

        System.out.println(v3.toString());
        System.out.println(v4.toString());

        try {
            v3.add(v4);
            fail("Vector (0,0,0) not valid");
        }
        catch  (IllegalArgumentException e)
        {
            assertTrue(e.getMessage()!= null);
        }
//        assertTrue(v3.length() >84);
        assertEquals(84,v3.length(),0.659);
    }

    @Test
    void length() {
    }

    @Test
    void testnormalize() {
        Vector v = new Vector(3.5, -5, 10);
        v.normalize();
        assertEquals(1, v.length(), 1e-10);

        try {
            Vector v1 = new Vector(0, 0, 0);
            v.normalize();
            fail("Didn't throw divide by zero exception!");
        } catch (IllegalArgumentException ex) {
            assertEquals("Point3D(0.0,0.0,0.0) not valid for vector head", ex.getMessage());
        }
        assertTrue(true);
    }
}
//
///**
// *
// */
//package unittests;
//
//        import static org.junit.Assert.*;
//        import org.junit.Test;
//
//        import static primitives.Util.*;
//        import primitives.*;
//
///**
// * Unit tests for primitives.Vector class
// *
// * @author Dan
// */
//public class VectorTests {
//    /**
//     * Test method for {@link primitives.Vector#add(primitives.Vector)}.
//     */
//    @Test
//    public void testAdd() {
//        // ============ Equivalence Partitions Tests ==============
//        // TC01: Simple test
//        assertEquals("Wrong vector add", new Vector(1, 1, 1),
//                new Vector(2, 3, 4).add(new Vector(-1, -2, -3)));
//
//        // =============== Boundary Values Tests ==================
//        // TC11: test adding v + (-v)
//        try {
//            new Vector(1, 2, 3).add(new Vector(-1, -2, -3));
//            fail("Add v plus -v must throw exception");
//        } catch (IllegalArgumentException e) {}
//    }
//
//    /**
//     * Test method for {@link primitives.Vector#subtract(primitives.Vector)}.
//     */
//    @Test
//    public void testSubtract() {
//        // ============ Equivalence Partitions Tests ==============
//        // TC01: Simple test
//        assertEquals("Wrong vector subtract", new Vector(1, 1, 1),
//                new Vector(2, 3, 4).subtract(new Vector(1, 2, 3)));
//
//        // =============== Boundary Values Tests ==================
//        // TC11: test subtracting same vector
//        try {
//            new Vector(1, 2, 3).subtract(new Vector(1, 2, 3));
//            fail("Subtract v from v must throw exception");
//        } catch (IllegalArgumentException e) {}
//    }
//
//    /**
//     * Test method for {@link primitives.Point3D#subtract(primitives.Point3D)}.
//     */
//    @Test
//    public void testPointSubtract() {
//        // ============ Equivalence Partitions Tests ==============
//        // TC01: Simple test
//        assertEquals("Wrong point subtract", new Vector(1, 1, 1),
//                new Point3D(2, 3, 4).subtract(new Point3D(1, 2, 3)));
//
//        // =============== Boundary Values Tests ==================
//        // TC11: test subtracting same point
//        try {
//            new Point3D(1, 2, 3).subtract(new Point3D(1, 2, 3));
//            fail("Subtract P from P must throw exception");
//        } catch (IllegalArgumentException e) {}
//    }
//
//    /**
//     * Test method for {@link primitives.Vector#scale(double)}.
//     */
//    @Test
//    public void testScale() {
//        // ============ Equivalence Partitions Tests ==============
//        // TC01: Simple test
//        assertEquals("Wrong vector scale", new Vector(2, 4, 6),
//                new Vector(1, 2, 3).scale(2));
//
//        // =============== Boundary Values Tests ==================
//        // TC11: test adding v + (-v)
//        try {
//            new Vector(1, 2, 3).scale(0d);
//            fail("Scale by 0 must throw exception");
//        } catch (IllegalArgumentException e) {}
//    }
//
//    /**
//     * Test method for {@link primitives.Vector#dotProduct(primitives.Vector)}.
//     */
//    @Test
//    public void testDotProduct() {
//        Vector v1 = new Vector(1, 2, 3);
//
//        // ============ Equivalence Partitions Tests ==============
//        // TC01: Simple dotProduct test
//        Vector v2 = new Vector(-2, -4, -6);
//        assertEquals("dotProduct() wrong value", -28d, v1.dotProduct(v2), 0.00001);
//
//        // =============== Boundary Values Tests ==================
//        // TC11: dotProduct for orthogonal vectors
//        Vector v3 = new Vector(0, 3, -2);
//        assertEquals("dotProduct() for orthogonal vectors is not zero", 0d, v1.dotProduct(v3), 0.00001);
//    }
//
//    /**
//     * Test method for {@link primitives.Vector#crossProduct(primitives.Vector)}.
//     */
//    @Test
//    public void testCrossProduct() {
//        Vector v1 = new Vector(1, 2, 3);
//
//        // ============ Equivalence Partitions Tests ==============
//        Vector v2 = new Vector(0, 3, -2);
//        Vector vr = v1.crossProduct(v2);
//
//        // TC01: Test that length of cross-product is proper (orthogonal vectors taken
//        // for simplicity)
//        assertEquals("crossProduct() wrong result length", v1.length() * v2.length(), vr.length(), 0.00001);
//
//        // TC02: Test cross-product result orthogonality to its operands
//        assertTrue("crossProduct() result is not orthogonal to 1st operand", isZero(vr.dotProduct(v1)));
//        assertTrue("crossProduct() result is not orthogonal to 2nd operand", isZero(vr.dotProduct(v2)));
//
//        // =============== Boundary Values Tests ==================
//        // TC11: test zero vector from cross-productof co-lined vectors
//        Vector v3 = new Vector(-2, -4, -6);
//        try {
//            v1.crossProduct(v3);
//            fail("crossProduct() for parallel vectors does not throw an exception");
//        } catch (Exception e) {}
//    }
//
//    /**
//     * Test method for {@link primitives.Vector#lengthSquared()}.
//     */
//    @Test
//    public void testLengthSquared() {
//        // ============ Equivalence Partitions Tests ==============
//        // TC01: Simple test
//        assertEquals("lengthSquared() wrong value", 14d, new Vector(1, 2, 3).lengthSquared(), 0.00001);
//    }
//
//    /**
//     * Test method for {@link primitives.Vector#length()}.
//     */
//    @Test
//    public void testLength() {
//        // TC01: Simple test
//        assertEquals("length() wrong value", 5d, new Vector(0, 3, 4).length(), 0.00001);
//    }
//
//    /**
//     * Test method for {@link primitives.Vector#normalized()}.
//     */
//    @Test
//    public void testNormalized() {
//        Vector v = new Vector(0, 3, 4);
//        Vector n = v.normalized();
//        // ============ Equivalence Partitions Tests ==============
//        // TC01: Simple test
//        assertFalse("normalized() changes the vector itself", v == n);
//        assertEquals("wrong normalized vector length", 1d, n.lengthSquared(), 0.00001);
//        try {
//            v.crossProduct(n);
//            fail("normalized vector is not in the same direction");
//        } catch (IllegalArgumentException e) {}
//        assertEquals("wrong normalized vector", new Vector(0, 0.6, 0.8), n);
//    }
//
//    /**
//     * Test method for {@link primitives.Vector#normalize()}.
//     */
//    @Test
//    public void testNormalize() {
//        Vector v = new Vector(0, 3, 4);
//        Vector n = v.normalize();
//        // ============ Equivalence Partitions Tests ==============
//        // TC01: Simple test
//        assertTrue("normalize() does not change the vector itself", v == n);
//        assertEquals("wrong normalized vector length", 1d, v.lengthSquared(), 0.00001);
//        assertEquals("wrong normalized vector", new Vector(0, 0.6, 0.8), v);
//    }
//
//}
