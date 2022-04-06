Usage example :

Poses of type Pose3, Ranges of type double

> Values values;

> values.insert(1, pose1);

> values.insert(2, pose2);

> values.insert(3, pose3);

> values.insert(4, pose4);

> LandmarkInitialization3D f;

> f.addRange(1, r1);

> f.addRange(2, r2);

> f.addRange(3, r3);

> f.addRange(4, r4);

> Point3 point = f.initialize(values);

For further help on usage, refer to the corresponding unit-test script. 
