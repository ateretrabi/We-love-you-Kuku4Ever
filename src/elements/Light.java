package elements;

import primitives.Color;
import primitives.Point3D;

public abstract class Light {
    /**
     * _intensity value, set to protected
     */
    protected Color _intensity;

    public Color getIntensity() {
        return new Color(_intensity);
    }
}
