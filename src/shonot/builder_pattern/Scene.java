package shonot.builder_pattern;

import elements.*;
import geometries.*;
import primitives.Color;


import java.util.LinkedList;
import java.util.List;

/**
 * using Builder and Fluent Interface
 */
public class Scene {
    private final String _name;
    private final Geometries _geometries;

    private Color _background;
    private Camera _camera;
    private double _distance;
    private AmbientLight _ambientLight;
    private List<LightSource> _lights = null;


    public AmbientLight getAmbientLight() {
        return _ambientLight;
    }

    public Camera getCamera() {
        return _camera;
    }

    public Geometries getGeometries() {
        return _geometries;
    }

    public double getDistance() {
        return _distance;
    }

    public Scene(SceneBuilder builder) {
        this._name = builder.name;
        this._geometries = new Geometries();
        this._camera = builder.camera;
        this._distance = builder.distance;
        this._background = builder.background;
        this._ambientLight = builder.ambientLight;
        builder.validateUserObject(this);
    }

    public Color getBackground() {
        return this._background;
    }


    public List<LightSource> getLightSources() {
        return _lights;
    }

    public void addGeometries(Intersectable... intersectables) {
        for (Intersectable i : intersectables) {
            _geometries.add(i);
        }
    }

    public void removeGeometries(Intersectable... intersectables) {
        for (Intersectable i : intersectables) {
            _geometries.remove(i);
        }
    }

    public void addLights(LightSource light) {
        if (_lights == null) {
            _lights = new LinkedList<>();
        }
        _lights.add(light);
    }

    public static class SceneBuilder {
        private final String name;
        private Color background;
        private Camera camera;
        private double distance;
        private AmbientLight ambientLight;

        public SceneBuilder(String name) {
            this.name = name;
        }

        public Scene.SceneBuilder addBackground(Color background) {
            this.background = background;
            return this;
        }

        public Scene.SceneBuilder addCamera(Camera camera) {
            this.camera = camera;
            return this;
        }

        public Scene.SceneBuilder addDistance(double distance) {
            this.distance = distance;
            return this;
        }

        public Scene.SceneBuilder addAmbientLight(AmbientLight ambientLight) {
            this.ambientLight = ambientLight;
            return this;
        }

        public Scene build() {
            Scene scene = new Scene(this);
            return scene;
        }

        private void validateUserObject(Scene scene) {
            //Do some basic validations to check
            //if user object does not break any assumption of system
        }
    }
}
