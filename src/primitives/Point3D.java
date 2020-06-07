package primitives;

/**
 * Class Point3D is the basic class representing a point in a
 * 3D system.
 *
 * @author Corona isDead
 */
public class Point3D {
    /**
     * _x coordinate on the X axis
     * _y coordinate on the Y axis
     * _z coordinate on the Z axis
     */
    public final static Point3D ZERO = new Point3D(0.0, 0.0, 0.0);
    Coordinate _x;
    Coordinate _y;
    Coordinate _z;

    public Point3D(Coordinate _x, Coordinate _y, Coordinate _z) {
        this._x = _x;
        this._y = _y;
        this._z = _z;
    }

    public Point3D(Point3D point3D) {
        this._x = new Coordinate(point3D._x);
        this._y = new Coordinate(point3D._y);
        this._z = new Coordinate(point3D._z);
    }


    public Point3D(double x, double y, double z) {
        this(new Coordinate(x), new Coordinate(y), new Coordinate(z));
    }

    public Coordinate getX() {
        return new Coordinate(_x);
    }

    public Coordinate getY() {
        return new Coordinate(_y);
    }

    public Coordinate getZ() {
        return new Coordinate(_z);
    }

    public double distanceSquared(Point3D otherPoint3D) {
        return ((otherPoint3D._x._coord - this._x._coord) * (otherPoint3D._x._coord - this._x._coord) +
                (otherPoint3D._y._coord - this._y._coord) * (otherPoint3D._y._coord - this._y._coord) +
                (otherPoint3D._z._coord - this._z._coord) * (otherPoint3D._z._coord - this._z._coord));
    }

    public double distance(Point3D otherPoint3D) {
        return Math.sqrt(distanceSquared(otherPoint3D));
    }

    public Point3D add(Vector vector) {
        return new Point3D(
                this._x._coord + vector._head._x._coord,
                this._y._coord + vector._head._y._coord,
                this._z._coord + vector._head._z._coord);
    }

    public Point3D subtract(Vector vector) {
        return new Point3D(
                this._x._coord - vector._head._x._coord,
                this._y._coord - vector._head._y._coord,
                this._z._coord - vector._head._z._coord);
    }

    public Vector subtract(Point3D otherPoint3D) {
        return new Vector(
                new Point3D(
                        this._x._coord - otherPoint3D._x._coord,
                        this._y._coord - otherPoint3D._y._coord,
                        this._z._coord - otherPoint3D._z._coord));
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Point3D point3D = (Point3D) o;
        return _x.equals(point3D._x) &&
                _y.equals(point3D._y) &&
                _z.equals(point3D._z);
    }

    @Override
    public String toString() {
        return "(" +
                _x +
                ", " + _y +
                ", " + _z +
                ')';
    }


}
