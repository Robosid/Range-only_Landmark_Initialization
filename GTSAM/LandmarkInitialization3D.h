
/**
 * @file  LandmarkInitialization3D.h
 * @brief Range-only 3D Initialization
 * Author : Sid Mahapatra @ Flanders Make
 * Status: Trilaterate, MDS, Bancroft with algo selection succesful | Cleaning pending
 */

#pragma once

#include <gtsam_unstable/dllexport.h>
#include <gtsam/nonlinear/NonlinearFactor.h>
#include <gtsam/inference/Key.h>
#include <gtsam/geometry/Pose3.h>

#include <list>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Core>
//#include <Eigen/SVD>
//#include "ceres/ceres.h"

/*
using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;
*/
namespace gtsam {

class LandmarkInitialization3D: public NoiseModelFactor {
 protected:
  struct Sphere3 {
    Sphere3(const Point3& p, double r) :
        center(p), radius(r) {
    }
    Point3 center;
    double radius;
  };
/*
  struct CostFunctor {

    // Initializing the constant parameters.

     CostFunctor(std::vector<double>& a_n, double d_n)
     : a_n_(a_n), d_n_(d_n) {}


      // Residual
      // \f$f_n(x) = ((a_{n_1} - x_1)^2 + (a_{n_2} - x_2)^2 + (a_{n_3} - x_3)^2 - d_n^2)\f$

     template <typename T>
     bool operator()(const T* const x,
                     T* residual) const {
       residual[0] = ((T(a_n_[0]) - x[0]) * (T(a_n_[0]) - x[0])
                    + (T(a_n_[1]) - x[1]) * (T(a_n_[1]) - x[1])
                    + (T(a_n_[2]) - x[2]) * (T(a_n_[2]) - x[2])
                    - T(d_n_) * T(d_n_));
       return true;
     }

       const std::vector<double>& a_n_;
       const double d_n_;
  };
*/
  typedef LandmarkInitialization3D This;

  std::vector<double> measurements_;  ///< Range measurements
  double variance_;  ///< variance on noise

 public:
  /** Default constructor: don't use directly */
  LandmarkInitialization3D() {
  }

  /**
   * Constructor
   * @param s standard deviation of range measurement noise
   */
  explicit LandmarkInitialization3D(double s) :
      NoiseModelFactor(noiseModel::Isotropic::Sigma(0, s)), variance_(s * s) {
  }

  ~LandmarkInitialization3D() override {
  }

  /// Add a range measurement to a pose with given key.
  void addRange(Key key, double measuredRange) {
    if(std::find(keys_.begin(), keys_.end(), key) != keys_.end()) {
      throw std::invalid_argument(
          "LandmarkInitialization3D::addRange: adding duplicate measurement for key.");
    }
    keys_.push_back(key);
    measurements_.push_back(measuredRange);
    size_t n = keys_.size();
    // Since we add the errors, the noise variance adds
    noiseModel_ = noiseModel::Isotropic::Variance(1, n * variance_);
    //std::cout<<variance_<<std::endl;
  }

  // Testable

  /** print */
  void print(const std::string& s = "",
      const KeyFormatter& keyFormatter = DefaultKeyFormatter) const override {
    std::cout << s << "LandmarkInitialization3D with " << size() << " measurements\n";
    NoiseModelFactor::print(s);
  }

  /** Check if two factors are equal */
  bool equals(const NonlinearFactor& f, double tol = 1e-9) const override {
    return false;
  }

  double dist(Point3 point1, Point3 point2) const{
    double myDistance=sqrt(pow(point1(0)-point2(0),2.0)+pow(point1(1)-point2(1),2.0)+pow(point1(2)-point2(2),2.0));
    return myDistance;
  }

  // factor interface

  /**
   * 3 methods to locate a point from at least four pose-range pairs
   * Raise runtime_error if not well defined.
   */

// Method 0: Using Google's Ceres Solver
/*
  Point3 trilaterate_ceres(const Values& x) const {
    // create n Spheres corresponding to measured range around each pose
    std::list<Sphere3> Spheres;
    Point3 location;
    size_t n = size();
    for (size_t j = 0; j < n; j++) {
      const Pose3& pose = x.at<Pose3>(keys_[j]);
      Spheres.push_back(Sphere3(pose.translation(), measurements_[j]));
    }
    // The variables to solve for
    std::vector<double> loc = {0.0, 0.0, 0.0};
    // Build the problem
    ceres::Problem problem;
    // Add residual term to the problem usin the autodiff
    for (std::list<Sphere3>::iterator i = Spheres.begin(); i != Spheres.end(); i++) {
      std::vector<double> pos = {i->center(0), i->center(1), i->center(2)};
      double dist = i->radius;
      ceres::CostFunction *cost_function =
      new ceres::AutoDiffCostFunction<CostFunctor, 1, 3>(
        new CostFunctor(pos, dist));
      problem.AddResidualBlock(cost_function,
                       NULL,
                       loc.data());
    }

      // Run the solver
      ceres::Solver::Options options;
      options.minimizer_progress_to_stdout = false;
      ceres::Solver::Summary summary;
      ceres::Solve(options, &problem, &summary);
      for (int i = 0; i < 3; i++) {
        location(i) =  loc[i];
      }

      return location;
}
*/

  Eigen::MatrixXd calculate_distance_matrix(Eigen::MatrixXd pose_positions, bool squared = true) const {

    Eigen::MatrixXd oper = (pose_positions.array() * pose_positions.array()).matrix();
    Eigen::MatrixXd oper1 = oper.rowwise().sum();
    Eigen::MatrixXd Adot = oper1.replicate(1, pose_positions.rows());
    Eigen::MatrixXd Bdot = Adot.transpose();
    Eigen::MatrixXd ABdot = pose_positions*(pose_positions.transpose());
    Eigen::MatrixXd D_squared =  (Adot + Bdot - (ABdot * 2));

    if(squared == true){
      return D_squared;
    }
    else{
      return D_squared.cwiseSqrt();
    }
  }

  void Kabsch(Eigen::MatrixXd Ps, Eigen::MatrixXd Qs, Eigen::MatrixXd *Rs, Eigen::MatrixXd *t)const {

    if ((Ps.rows() != Qs.rows()) || (Ps.cols() != Qs.cols()))
    {
      std::cout<<"Ps and Qs shapes mismatch";
      throw std::runtime_error("Logic error");
    }
    Eigen::MatrixXd p = Ps.colwise().mean();
    Eigen::MatrixXd q = Qs.colwise().mean();
    Eigen::MatrixXd p0 = p.replicate(Ps.rows() ,1);
    Eigen::MatrixXd q0 = q.replicate(Qs.rows() ,1);
    Eigen::MatrixXd P_centered = Ps - p0;
    Eigen::MatrixXd Q_centered = Qs - q0;
    Eigen::MatrixXd Cs = (P_centered.transpose()) * Q_centered;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cs, Eigen::ComputeThinU | Eigen::ComputeThinV);
    *Rs = svd.matrixU() * (svd.matrixV()).transpose();
    *t = q - (p * (*Rs));
  }

  // Method 1: Using Multi-Dimensional Scaling
  Point3 MDS_loc_engine_3D(const Values& x)const  {

    Point3 location;
    std::list<Sphere3> Spheres;
    size_t s = size();
    if (keys_.size() != measurements_.size())
    {
      std::cout<<"Pose keys size and ranges size mismatch";
      throw std::runtime_error("Logic error");
    }

    for (size_t j = 0; j < s; j++) {
      const Pose3& pose = x.at<Pose3>(keys_[j]);
      Spheres.push_back(Sphere3(pose.translation(), measurements_[j]));
    }

    Eigen::MatrixXd pose_positions(s, 3);
    int count = 0;
    for (std::list<Sphere3>::const_iterator i = Spheres.begin(); i != Spheres.end(); ++i) {
      for(int j = 0; j < 3; ++j){         // for 3D
      pose_positions(count,j) = i->center(j);}
      count++;
    }

    Eigen::MatrixXd ranges(s, 1);
    count = 0;
    for (std::list<Sphere3>::const_iterator i = Spheres.begin(); i != Spheres.end(); ++i) {
      ranges(count, 0) = i->radius;
      count++;
    }

    int n = pose_positions.rows(), p = pose_positions.cols(), m = ranges.cols();

    Eigen::MatrixXd D_squared = calculate_distance_matrix(pose_positions);
    Eigen::MatrixXd D_full = Eigen::MatrixXd::Zero(D_squared.rows()+1, D_squared.cols()+1);
    D_full.block(0,0 , D_squared.rows(), D_squared.cols()) = D_squared.block(0,0 , D_squared.rows(), D_squared.cols());
    D_full.block(0, D_squared.cols(), D_squared.rows(), 1) = ranges.cwiseAbs2();
    D_full.block(D_squared.rows(), 0, 1, ranges.rows()) = (ranges.transpose()).cwiseAbs2();
    Eigen::MatrixXd Cs = Eigen::MatrixXd::Identity(n+m, n+m) - Eigen::MatrixXd::Constant(n+m,n+m, 1.0)/(n+m);
    Eigen::MatrixXd Bs = (-0.5 * Cs) * D_full * Cs;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Bs, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd Sx = (svd.singularValues()).block(0,0,p,1);
    Eigen::MatrixXd Dia = Eigen::MatrixXd::Zero(p,p);
    for (int i = 0; i < p; i++){
      for (int j = 0; j < p; j++){
        if (i == j)
          Dia(i,j) = std::sqrt(Sx(i));
      }
    }
    Eigen::MatrixXd Ps = ((svd.matrixU()).block(0,0,(svd.matrixU()).rows(), p)) * Dia;
    Eigen::MatrixXd Rs, t;
    Eigen::MatrixXd out;
    Kabsch(Ps.block(0, 0, Ps.rows()-1, Ps.cols()), pose_positions, &Rs, &t);
    out = ((Ps.block(Ps.rows()-1, 0, 1, Ps.cols())) * Rs) + t;

    for (int i = 0; i < out.cols(); i++)
      location(i) = out(i);

    return location;
  }

  // dotLorentz with dimension = 0
  Eigen::MatrixXd dotLorentz(Eigen::MatrixXd u, Eigen::MatrixXd v) const {

    Eigen::MatrixXd out;
    Eigen::MatrixXd helper;
    Eigen::MatrixXd assist;

    helper = (u.array() * v.array()).matrix();
    assist = (helper.block(0, 0, helper.rows(), helper.cols()-1)).rowwise().sum();
    out = assist - helper.block(0, helper.cols()-1, helper.rows(), 1);
    return out;
  }

  // method for calculating the pseudo-Inverse as recommended by Eigen developers
  template<typename _Matrix_Type_>
  _Matrix_Type_ pseudoInverse(const _Matrix_Type_ &a, double epsilon = std::numeric_limits<double>::epsilon()) const {
          //For a square matrix
  	//Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeFullU | Eigen::ComputeFullV);
          // For a non-square matrix
    Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeThinU | Eigen::ComputeThinV);
  	double tolerance = epsilon * std::max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);
  	return svd.matrixV() * (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
  }

  // Method 2: Using Bancroft Algorithm
  Point3 Bancroft_loc_engine_3D(const Values& x) const {

    Point3 location;
    std::list<Sphere3> Spheres;
    size_t s = size();
    if (keys_.size() != measurements_.size())
    {
      std::cout<<"Pose keys size and ranges size mismatch";
      throw std::runtime_error("Logic error");
    }

    for (size_t j = 0; j < s; j++) {
      const Pose3& pose = x.at<Pose3>(keys_[j]);
      Spheres.push_back(Sphere3(pose.translation(), measurements_[j]));
    }

    Eigen::MatrixXd pose_positions(s, 3);
    int count = 0;
    for (std::list<Sphere3>::const_iterator i = Spheres.begin(); i != Spheres.end(); ++i) {
      for(int j = 0; j < 3; ++j){         // for 3D
      pose_positions(count,j) = i->center(j);}
      count++;
    }

    Eigen::MatrixXd ranges(s, 1);
    count = 0;
    for (std::list<Sphere3>::const_iterator i = Spheres.begin(); i != Spheres.end(); ++i) {
      ranges(count, 0) = i->radius;
      count++;
    }

    Eigen::MatrixXd Bs = Eigen::MatrixXd::Zero(pose_positions.rows(), pose_positions.cols()+1);
    Bs.block(0, 0, pose_positions.rows(), pose_positions.cols()) = pose_positions;
    Bs.block(0, pose_positions.cols(), ranges.rows(), 1) = ranges * (-1);
    Eigen::MatrixXd Bs_ = pseudoInverse(Bs);
    Eigen::MatrixXd a = (dotLorentz(Bs, Bs) * 0.5);
    double c2 = (dotLorentz((Bs_.rowwise().sum()).transpose(), (Bs_.rowwise().sum()).transpose()))(0);
    double c1 = ((dotLorentz((Bs_ * a), (Bs_.rowwise().sum()).transpose()))(0) - 1) * 2;
    double c0 = (dotLorentz((Bs_ * a).transpose(), (Bs_ * a).transpose()))(0);
    double root1, root2;
    if ((c1 * c2 - 4 * c0 * c2) >= 0){
      root1 = (-c1 + std::sqrt(c1*c1-4*c0*c2)) / (2*c2);
      root2 = (-c1 - std::sqrt(c1*c1-4*c0*c2)) / (2*c2);
    }
    else{
      root1 = -c1 / (2*c2);
      root2 = root1;
    }
    Eigen::MatrixXd u1 = (Bs_ * (a + (Eigen::MatrixXd::Ones(a.rows(), a.cols()) * root1))).transpose();
    Eigen::MatrixXd u2 = (Bs_ * (a + (Eigen::MatrixXd::Ones(a.rows(), a.cols()) * root2))).transpose();;
    Eigen::MatrixXd Us(u1.rows() + u2.rows(), u1.cols());
    Us << u1,u2;
    Us = Us.cwiseAbs();
    if (Us(0, u1.cols()-1) < Us(1, u2.cols()-1)){
      for (int i=0; i < u1.cols() - 1; i++)
        location(i) = u1(0, i);
    }
    else{
      for (int i=0; i < u2.cols() - 1; i++)
        location(i) = u2(0, i);
    }

    return location;
  }

  // Method 3: Using Trilateration with SVD
  Point3 trilaterate(const Values& x) const {
    // create n Spheres corresponding to measured range around each pose
    std::list<Sphere3> Spheres;
    size_t n = size();
    Point3 location;
    for (size_t j = 0; j < n; j++) {
      const Pose3& pose = x.at<Pose3>(keys_[j]);
      Spheres.push_back(Sphere3(pose.translation(), measurements_[j]));
    }

    // Define the matrix that we are going to use
    size_t rows = n * (n - 1) / 2;
    Eigen::MatrixXd m(rows, 3);
    Eigen::VectorXd b(rows);

    // Fill in matrix according to the equations
    size_t row = 0;
    double x1, x2, y1, y2, z1, z2, r1, r2;
    for (std::list<Sphere3>::const_iterator i = Spheres.begin(); i != Spheres.end(); ++i) {
      for (std::list<Sphere3>::const_iterator j = i; j != Spheres.end(); ++j) {
        // distance between Sphere centers
        double d = distance3(i->center, j->center);
        // skip Spheres that are in the same location
        if (d < 1e-9)
          continue;
        x1 = i->center(0), y1 = i->center(1), z1 = i->center(2);
        x2 = j->center(0), y2 = j->center(1), z2 = j->center(2);
        r1 = i->radius;
        r2 = j->radius;
        m(row, 0) = x1 - x2;
        m(row, 1) = y1 - y2;
        m(row, 2) = z1 - z2;
        b(row) = ((pow(x1, 2)-pow(x2, 2)) +
                  (pow(y1, 2)-pow(y2, 2)) +
                  (pow(z1, 2)-pow(z2, 2)) -
                  (pow(r1, 2) - pow(r2, 2))) / 2;
        row++;
      }
    }

    // Then calculate to solve the equations, using the least square solution
    try {
      location = m.jacobiSvd(Eigen::ComputeThinU|Eigen::ComputeThinV).solve(b);
      return location;
  }
    catch(...){
      std::cout<<"Error: Trilateration failed";
      throw std::runtime_error("Runtime error");
  }
    }

  bool checkInitialization(const Values& x, Point3 loc_estimate, float pdopths = 20, float sigmaths = 0.1, float planeratioths = 1) const{

      std::list<Sphere3> Spheres;
      size_t s = size();
      if (keys_.size() != measurements_.size())
      {
        std::cout<<"Pose keys size and ranges size mismatch";
        throw std::runtime_error("Logic error");
      }

      for (size_t j = 0; j < s; j++) {
        const Pose3& pose = x.at<Pose3>(keys_[j]);
        Spheres.push_back(Sphere3(pose.translation(), measurements_[j]));
      }

      Eigen::MatrixXd pose_positions(s, 3);
      int count = 0;
      for (std::list<Sphere3>::const_iterator i = Spheres.begin(); i != Spheres.end(); ++i) {
        for(int j = 0; j < 3; ++j){         // for 3D
        pose_positions(count,j) = i->center(j);}
        count++;
      }

      Eigen::MatrixXd ranges(s, 1);
      count = 0;
      for (std::list<Sphere3>::const_iterator i = Spheres.begin(); i != Spheres.end(); ++i) {
        ranges(count, 0) = i->radius;
        count++;
      }

      Eigen::MatrixXd pos_est(1, 3);
      for (int i = 0; i < 3; i++)
        pos_est(0, i) = loc_estimate(i);

      int n = pose_positions.rows();
      pdopths = pdopths / sqrt(n);
      Eigen::MatrixXd As = pose_positions - pos_est.replicate(pose_positions.rows(), 1);
      for (int i = 0; i < As.rows(); i++){
          for (int j = 0; j < As.cols(); j++){
            As(i,j) = As(i,j) / ranges(i, 0);
          }
      }

      Eigen::MatrixXd Qs = (As.transpose() * As).inverse();
      float pdop = sqrt(Qs.diagonal().sum());
      Eigen::MatrixXd ranges_estimate = (pose_positions - pos_est.replicate(pose_positions.rows(), 1)).rowwise().norm();
      float sigma = (ranges_estimate - ranges).norm() / sqrt(n-1);
      Eigen::MatrixXd mean_poses = pose_positions.colwise().mean();
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(pose_positions - mean_poses.replicate(pose_positions.rows(),1), Eigen::ComputeThinU | Eigen::ComputeThinV);
      Eigen::MatrixXd Vs = svd.matrixV();
      for(int i = 0; i<Vs.rows(); i++)
        Vs(i,0) = Vs(i,0) * -1;
      Eigen::MatrixXd q = pos_est - mean_poses;
      double angle = asin((abs((Vs.block(0, Vs.cols()-1, Vs.rows(), 1).adjoint() * q.transpose())(0,0))) / q.norm());
      double b = -q.norm() * sin(angle) + sqrt(pow(sin(angle),2) * pow(q.norm(),2) + sigma * (2 * q.norm() + sigma));
      double plane_ratio = (svd.singularValues()(svd.singularValues().rows()-1, 0)) / sqrt(n) / b;
      bool out = (plane_ratio > planeratioths) && (pdop < pdopths) && (sigma < sigmaths);

      return out;

  }

  Point3 initialize(const Values& x) const {

    Point3 optimizedPoint = MDS_loc_engine_3D(x);
    if (checkInitialization(x, optimizedPoint) == true)
      return optimizedPoint;
    else {
      optimizedPoint = Bancroft_loc_engine_3D(x);
      if (checkInitialization(x, optimizedPoint) == true)
        return optimizedPoint;
      else {
        optimizedPoint = trilaterate(x);
        if (checkInitialization(x, optimizedPoint) == true)
          return optimizedPoint;
        else
          return Point3 (-1, -1, -1);
        }
     }
  }

  // FOR DEBUG ONLY
  Vector unwhitenedError(const Values& x,
      boost::optional<std::vector<Matrix>&> H = boost::none) const override {
    size_t n = size();
    if (n < 4) {
      if (H) {
        // set Jacobians to zero for n<4
        // To locate position in a 3D space, have to get at least 4 becaons
        for (size_t j = 0; j < n; j++)
          (*H)[j] = Matrix::Zero(3, 6); // Default: 3-6, TRY n ???
      }
      return Z_1x1;
    } else {
      Vector error = Z_1x1;

      // TODO: Should we have a (better?) variant that does this in relative coordinates ?
      Point3 optimizedPoint = MDS_loc_engine_3D(x);
      bool check = checkInitialization(x, optimizedPoint);
      std::cout<<check;
    //  Point3 optimizedPoint = trilaterate(x);

      // TODO: triangulation should be followed by an optimization given poses
      // now evaluate the errors between predicted and measured range
      for (size_t j = 0; j < n; j++) {
        const Pose3& pose = x.at<Pose3>(keys_[j]);
        if (H)
          // also calculate 1*6 derivative for each of the n poses
          error[0] += pose.range(optimizedPoint, (*H)[j]) - measurements_[j];
        else
          error[0] += dist(pose.translation(), optimizedPoint) - measurements_[j];
      }
      return error;
   }
}

  /// @return a deep copy of this factor
  gtsam::NonlinearFactor::shared_ptr clone() const override {
    return boost::static_pointer_cast<gtsam::NonlinearFactor>(
        gtsam::NonlinearFactor::shared_ptr(new This(*this)));
  }
};
}  // \namespace gtsam
