
/**
 *  @file  testLandmarkInitialization3D.cpp
 *  @brief Unit tests for LandmarkInitialization3D Class
 */

#include <gtsam/sam/LandmarkInitialization3D.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <CppUnitLite/TestHarness.h>

using namespace std;
using namespace gtsam;

static const double sigma = 2.0;

// Test situation:
static const Point3 p(4, 5, 2);
static const Pose3 pose1(Rot3::Rodrigues(0,0,0), Point3(0,0,0)), pose2(Rot3::Rodrigues(0,0,0), Point3(2,0,0)), pose3(Rot3::Rodrigues(0,0,0), Point3(0,2,0)), pose4(Rot3::Rodrigues(0,0,0), Point3(1,1,2));
//static const double r1 = pose1.range(p), r2 = pose2.range(p), r3 = pose3.range(p), r4 = pose4.range(p);

/* ************************************************************************* */

TEST( LandmarkInitialization3D, constructor ) {
  LandmarkInitialization3D f1;
  LONGS_EQUAL(0, f1.size())
  LandmarkInitialization3D f2(sigma);
  LONGS_EQUAL(0, f2.size())
}

/* ************************************************************************* */

TEST( LandmarkInitialization3D, addRange ) {
  LandmarkInitialization3D f(sigma);
  f.addRange(1, 10);
  f.addRange(2, 12);
  LONGS_EQUAL(2, f.size())
}

/* ************************************************************************* */
/*
TEST( LandmarkInitialization3D, scenario ) {
  //DOUBLES_EQUAL(expected, actual, tolerance)
  DOUBLES_EQUAL(13.34, r1, 1e-9);
  DOUBLES_EQUAL(sqrt(100.0+25.0), r2, 1e-9);
  DOUBLES_EQUAL(sqrt(50.0), r3, 1e-9);
}
*/

/* ************************************************************************* */

TEST( LandmarkInitialization3D, unwhitenedError ) {
  Values values; // all correct
  values.insert(1, pose1);
  values.insert(2, pose2);
  values.insert(3, pose3);
  values.insert(4, pose4);

  LandmarkInitialization3D f(99.9);
  double r1 = f.dist(pose1.translation(), p), r2 = f.dist(pose2.translation(), p), r3 = f.dist(pose3.translation(), p), r4 = f.dist(pose4.translation(), p);
  f.addRange(1, r1);

  // Check Jacobian for n==1
  vector<Matrix> H1(1);
  f.unwhitenedError(values, H1); // with H now !
  CHECK(assert_equal(Matrix::Zero(3,6), H1.front()));

  // Whenever there are three ranges or less, error should be zero
  Vector actual1 = f.unwhitenedError(values);
  EXPECT(assert_equal((Vector(1) << 0.0).finished(), actual1));
  f.addRange(2, r2);
  Vector actual2 = f.unwhitenedError(values);
  EXPECT(assert_equal((Vector(1) << 0.0).finished(), actual2));
  f.addRange(3, r3);
  Vector actual3 = f.unwhitenedError(values);
  EXPECT(assert_equal((Vector(1) << 0.0).finished(), actual3));

  //Point3 check_point0 = f.trilaterate(values);
  //CHECK(assert_equal(p, check_point0));

  f.addRange(4, r4);
  Vector actual4 = f.unwhitenedError(values);
  //EXPECT_LONGS_EQUAL(4, f.keys().size());
  //EXPECT(assert_equal((Vector(1) << 0.0).finished(), actual4));

  //f.addRange(5, r5);
  //EXPECT_LONGS_EQUAL(5, f.keys().size());
  //Point3 check_point = f.MDS_loc_engine_3D(values);
  Point3 check_point = f.Bancroft_loc_engine_3D(values);
  CHECK(assert_equal(p, check_point));

  // Check keys and Jacobian
  //vector<Matrix> H(4);
  //Vector actual5 = f.unwhitenedError(values, H); // with H now !
  //CHECK(assert_equal((Vector(1) << 0.0).finished(), actual5));
  //CHECK(assert_equal((Matrix(1, 6) << 0.0,-1.0,0.0).finished(), H.front()));
  //CHECK(assert_equal((Matrix(1, 6) << sqrt(2.0)/2,-sqrt(2.0)/2,0.0).finished(), H.back()));

  // Test clone
  NonlinearFactor::shared_ptr clone = f.clone();
  EXPECT_LONGS_EQUAL(4, clone->keys().size());
}

/* ************************************************************************* */
/*
TEST( LandmarkInitialization3D, optimization ) {
  LandmarkInitialization3D f(sigma);
  f.addRange(1, r1);
  f.addRange(2, r2);
  f.addRange(3, r3);

  // Create initial value for optimization
  Values initial;
  initial.insert(1, pose1);
  initial.insert(2, pose2);
  initial.insert(3, Pose2(5, 6, 0)); // does not satisfy range measurement
  Vector actual5 = f.unwhitenedError(initial);
  EXPECT(assert_equal((Vector(1) << sqrt(25.0+16.0)-sqrt(50.0)).finished(), actual5));

  // Create Factor graph
  NonlinearFactorGraph graph;
  graph.push_back(f);
  const noiseModel::Base::shared_ptr //
  priorNoise = noiseModel::Diagonal::Sigmas(Vector3(1, 1, M_PI));
  graph.addPrior(1, pose1, priorNoise);
  graph.addPrior(2, pose2, priorNoise);

  // Try optimizing
  LevenbergMarquardtParams params;
  //  params.setVerbosity("ERROR");
  LevenbergMarquardtOptimizer optimizer(graph, initial, params);
  Values result = optimizer.optimize();
  EXPECT(assert_equal(initial.at<Pose2>(1), result.at<Pose2>(1)));
  EXPECT(assert_equal(initial.at<Pose2>(2), result.at<Pose2>(2)));
  // only the third pose will be changed, converges on following:
  EXPECT(assert_equal(Pose2(5.52159, 5.582727, 0), result.at<Pose2>(3),1e-5));
}
*/
/* ************************************************************************* */
int main() {
  TestResult tr;
  return TestRegistry::runAllTests(tr);
}
/* ************************************************************************* */
