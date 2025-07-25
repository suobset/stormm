// -*-c++-*-
#ifndef STORMM_WATCHER_H
#define STORMM_WATCHER_H

#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/behavior.h"
#include "DataTypes/stormm_vector_types.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"

namespace stormm {
namespace review {

using card::Hybrid;
using card::HybridTargetLevel;
using constants::ExceptionResponse;
using synthesis::AtomGraphSynthesis;
using synthesis::PhaseSpaceSynthesis;

/// \brief The maximum number of reports that a Watcher object will hold
constexpr int default_max_watcher_reports = 1000;

/// \brief The default threshold at which to report forces, in kcal/mol-A
constexpr int default_watcher_force_threshold = 128.0;

/// \brief The default particle speed at which to report high velocities, in kcal/mol-A
constexpr int default_watcher_speed_threshold = 1.0;

/// \brief The minimum force threshold for reporting legitimately suspicious forces (this is more
///        of a guard against bogus values, not bad user input)
constexpr float minimum_force_threshold = 1.0;

/// \brief The minimum particle speed threshold for reporting legitimately suspicious forces (this
///        is more of a guard against bogus values, not bad user input)
constexpr float minimum_speed_threshold = 1.0;

/// \brief The writeable abstract for a Watcher class will be provided as a formal argument to
///        various molecular simulations functions that record checks for it.
struct WatcherWriter {
public:
  
  /// \brief As with other abstracts, the constructor takes a list of all relevant pointers and
  ///        critical constants.
  WatcherWriter(int nsystem_in, int max_reports_in, float force_limit_in, float speed_limit_in,
                bool track_purge_in, int* nforce_in, int* nspeed_in, int* nrattle_in,
                int* nshake_in, float4* forces_in, int* force_steps_in, int* force_stages_in,
                float4* speeds_in, int* speed_steps_in, int* speed_stages_in,
                uint2* rattle_fails_in, uint2* shake_fails_in, float* rattle_ext_in,
                float* shake_ext_in, float* xvel_purge_in, float* yvel_purge_in,
                float* zvel_purge_in, float* xang_purge_in, float* yang_purge_in,
                float* zang_purge_in);

  /// \brief Like most abstracts, the presence of const members makes copy and move assignment
  ///        impossible except in very recent C++ implementations.  The copy and move constructors
  ///        are legal and can be taken in their default forms.
  ///
  /// \param original  The original object to copy or move
  /// \{
  WatcherWriter(const WatcherWriter &original) = default;
  WatcherWriter(WatcherWriter &&original) = default;
  /// \}
  
  const int nsystem;        ///< The number of systems tracked by the Watcher object
  const int max_reports;    ///< The maximum number of reports to present in any given category
  const float force_limit;  ///< The largest force (magnitude) that will not be reported
  const float speed_limit;  ///< The largest particle speed (velocity magnitude) that will not be
                            ///<   reported
  const bool track_purge;   ///< Indicate whether the Watcher will also track momentum purges, up
                            ///<   to a certain number for each system.  This is using the Watcher
                            ///<   in a capacity rather like the ScoreCard for recording energetic
                            ///<   state variables.
  int* nforce;              ///< The number of forces observed to violate the stated bounds
  int* nspeed;              ///< The number of particles observed moving faster than has been
                            ///<   deemed acceptable
  int* nrattle;             ///< The number of unconverged RATTLE executions observed over the
                            ///<   course of the simulation
  int* nshake;              ///< The number of unconverged SHAKE executions observed over the
                            ///<   course of the simulation
  float4* forces;           ///< Cartesian components and atom indices of large forces
  int* force_steps;         ///< Steps on which each large force is observed
  int* force_stages;        ///< Stages of the integration cycle (this is not the WHITE / BLACK
                            ///<   coordinate cycle) in which each large force is observed
  float4* speeds;           ///< Cartesian components and atom indices of high speeds observed in
                            ///<   the simulation
  int* speed_steps;         ///< Steps on which each high speed is observed
  int* speed_stages;        ///< Stages of the integration cycle at which each high speed is
                            ///<   observed
  uint2* rattle_fails;      ///< Details of RATTLE group convergence failures
  uint2* shake_fails;       ///< Details of SHAKE group convergence failures
  float* rattle_ext;        ///< The extent to which each RATTLE group failed to converge.  How far
                            ///<   from the convergence tolerance did the iterations stop?
  float* shake_ext;         ///< The extent to which each SHAKE group failed to converge
  float* xvel_purge;        ///< A compendium of Cartesian X velocity purges for each system,
                            ///<   ordered as the first purge for systems 0, 1, ..., N (where N is
                            ///<   the number of systems), the second purge for systems 0, 1, ...,
                            ///<   N, up to the Mth purge for all systems
  float* yvel_purge;        ///< A compendium of Cartesian Y velocity purges for each system
  float* zvel_purge;        ///< A compendium of Cartesian Z velocity purges for each system
  float* xang_purge;        ///< Angular momentum purges about the X axis for each system
  float* yang_purge;        ///< Angular momentum purges about the Y axis for each system
  float* zang_purge;        ///< Angular momentum purges about the Z axis for each system
};

/// \brief The read-only abstract for the Watcher class is available to make reports out of a
///        Watcher object passed by const reference.
struct WatcherReader {
public:

  /// \brief As with other abstracts, the constructor takes a list of all relevant pointers and
  ///        critical constants.
  /// \{
  WatcherReader(int nsystem_in, int max_reports_in, float force_limit_in, float speed_limit_in,
                bool track_purge_in, const int* nforce_in, const int* nspeed_in,
                const int* nrattle_in, const int* nshake_in, const float4* forces_in,
                const int* force_steps_in, const int* force_stages_in, const float4* speeds_in,
                const int* speed_steps_in, const int* speed_stages_in,
                const uint2* rattle_fails_in, const uint2* shake_fails_in,
                const float* rattle_ext_in, const float* shake_ext_in, const float* xvel_purge_in,
                const float* yvel_purge_in, const float* zvel_purge_in, const float* xang_purge_in,
                const float* yang_purge_in, const float* zang_purge_in);

  WatcherReader(const WatcherWriter &w);
  /// \}

  /// \brief Like most abstracts, the presence of const members makes copy and move assignment
  ///        impossible except in very recent C++ implementations.  The copy and move constructors
  ///        are legal and can be taken in their default forms.
  ///
  /// \param original  The original object to copy or move
  /// \{
  WatcherReader(const WatcherReader &original) = default;
  WatcherReader(WatcherReader &&original) = default;
  /// \}

  const int nsystem;          ///< The number of systems tracked by the Watcher object
  const int max_reports;      ///< The maximum number of reports to present in any given category
  const float force_limit;    ///< The largest force (magnitude) that will not be reported
  const float speed_limit;    ///< The largest particle speed (velocity magnitude) that will not be
                              ///<   reported
  const bool track_purge;     ///< Indicate whether the Watcher will also track momentum purges, up
                              ///<   to a certain number for each system.  This is using the
                              ///<   Watcher in a capacity rather like the ScoreCard for recording
                              ///<   energetic state variables.
  const int* nforce;          ///< The number of forces observed to violate the stated bounds
  const int* nspeed;          ///< The number of particles observed moving faster than has been
                              ///<   deemed acceptable
  const int* nrattle;         ///< The number of unconverged RATTLE executions observed over the
                              ///<   course of the simulation
  const int* nshake;          ///< The number of unconverged SHAKE executions observed over the
                              ///<   course of the simulation
  const float4* forces;       ///< Cartesian components and atom indices of large forces
  const int* force_steps;     ///< Steps on which each large force is observed
  const int* force_stages;    ///< Stages of the integration cycle (this is not the WHITE / BLACK
                              ///<   coordinate cycle) in which each large force is observed
  const float4* speeds;       ///< Cartesian components and atom indices of high speeds observed in
                              ///<   the simulation
  const int* speed_steps;     ///< Steps on which each high speed is observed
  const int* speed_stages;    ///< Stages of the integration cycle at which each high speed is
                              ///<   observed
  const uint2* rattle_fails;  ///< Details of RATTLE group convergence failures
  const uint2* shake_fails;   ///< Details of SHAKE group convergence failures
  const float* rattle_ext;    ///< The extent to which each RATTLE group failed to converge.  How
                              ///<   far from the convergence tolerance did the iterations stop?
  const float* shake_ext;     ///< The extent to which each SHAKE group failed to converge
  const float* xvel_purge;    ///< A compendium of Cartesian X velocity purges for each system,
                              ///<   ordered as the first purge for systems 0, 1, ..., N (where N
                              ///<   is the number of systems), the second purge for systems 0, 1,
                              ///<   ..., N, up to the Mth purge for all systems
  const float* yvel_purge;    ///< A compendium of Cartesian Y velocity purges for each system
  const float* zvel_purge;    ///< A compendium of Cartesian Z velocity purges for each system
  const float* xang_purge;    ///< Angular momentum purges about the X axis for each system
  const float* yang_purge;    ///< Angular momentum purges about the Y axis for each system
  const float* zang_purge;    ///< Angular momentum purges about the Z axis for each system
};

/// \brief The Watcher class exists to take note of 'significant' events that may occur while
///        propagating a trajectory or during an energy minimization.  Various O(N) processes may
///        be recorded without adding significantly to the simulation cost, which could provide
///        valuable information as to why it may be failing or even indicate problems in situations
///        that appear otherwise normal.
class Watcher {
public:

  /// \brief The constructor requires a coordinate synthesis in order to allocate the proper amount
  ///        of space for the number of systems in the calulation.
  /// \{
  Watcher(const PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag,
          float force_threshold_in = default_watcher_force_threshold,
          float speed_threshold_in = default_watcher_speed_threshold,
          bool track_momentum_purge_in = false, int max_reports_in = default_max_watcher_reports,
          ExceptionResponse policy_in = ExceptionResponse::WARN);

  Watcher(const PhaseSpaceSynthesis &poly_ps, const AtomGraphSynthesis &poly_ag,
          float force_threshold_in = default_watcher_force_threshold,
          float speed_threshold_in = default_watcher_speed_threshold,
          bool track_momentum_purge_in = false, int max_reports_in = default_max_watcher_reports,
          ExceptionResponse policy_in = ExceptionResponse::WARN);
  /// \}

  /// \brief The presence of POINTER-kind Hybrid objects implies pointers to repair and thus
  ///        manual coding of copy and move constructors as well as copy and move assignment
  ///        operators.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object placed on the right-hand-side of the assignment operation
  /// \{
  Watcher(const Watcher &original);
  Watcher(Watcher &&original);
  Watcher& operator=(const Watcher &original);
  Watcher& operator=(Watcher &&original);
  /// \}

  /// \brief Get the maximum number of reports that the object is designed to hold.
  int getReportCount() const;
  
  /// \brief Get the number of systems from the associated coordinate synthesis.
  int getSystemCount() const;
  
  /// \brief Get the threshold at which the object will report large forces on atoms.
  float getForceThreshold() const;

  /// \brief Get the number of large forces on atoms found throughout the simulation.
  ///
  /// \param tier  Take the result from CPU host or GPU device memory
  int getLargeForceCount(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the list of large forces observed over the course of simulations on a set of
  ///        systems.
  ///
  /// \param tier  Take the result from CPU host or GPU device memory
  std::vector<float4> getLargeForces(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the list of step counts at which each large force in the list of all large forces
  ///        was observed.
  ///
  /// \param tier  Take the result from CPU host or GPU device memory
  std::vector<int> getLargeForceSteps(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the points in the integration cycle at which each of the large forces were
  ///        observed.  The integration cycle is not the coordinate time cycle, rather a more
  ///        detailed list of the stages within each time step that will help narrow down the
  ///        source of jarring forces or sudden movement.
  ///
  /// \param tier  Take the result from CPU host or GPU device memory
  std::vector<DynamicsStepStage>
  getLargeForceStages(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the threshold at which the object will report large velocities on atoms.
  float getSpeedThreshold() const;
  
  /// \brief Get the number of high velocities encountered in the simulation.
  ///
  /// \param tier  Take the result from CPU host or GPU device memory
  int getHighSpeedCount(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the list of high velocities observed over the course of simulations on a set of
  ///        systems.
  ///
  /// \param tier  Take the result from CPU host or GPU device memory
  std::vector<float4> getHighSpeeds(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the list of step counts at which each high velocity in the list of all high
  ///        veloicities was observed.
  ///
  /// \param tier  Take the result from CPU host or GPU device memory
  std::vector<int> getHighSpeedSteps(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the points in the integration cycle at which each high particle velocity was
  ///        found.
  ///
  /// \param tier  Take the result from CPU host or GPU device memory
  std::vector<DynamicsStepStage>
  getHighSpeedStages(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  
  /// \brief Get the number of times the RATTLE constraint iteration limit was reached (without the
  ///        requested convergence being achieved).
  ///
  /// \param tier          Obtain results from the CPU host or GPU device memory
  int getFailedRattleCount(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the number of times the SHAKE constraint iteration limit was reached (without the
  ///        requested convergence being achieved).  Descriptions of input parameters follow from
  ///        getRattleIterationViolations(), above.
  ///
  /// \param tier  Take the result from CPU host or GPU device memory
  int getFailedShakeCount(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  
  /// \brief Get the list of RATTLE violations--group IDs in the "x" member of the tuple (as listed
  ///        in the associated topology synthesis) and the step number at which the constraint
  ///        failed to converge in the "y" member.
  ///
  /// \param tier  Obtain results from the CPU host or GPU device memory
  std::vector<uint2> getRattleFailures(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the list of SHAKE violations--group IDs in the "x" member of the tuple (as listed
  ///        in the associated topology synthesis) and the step number at which the constraint
  ///        failed to converge in the "y" member.
  ///
  /// \param tier  Obtain results from the CPU host or GPU device memory
  std::vector<uint2> getShakeFailures(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the violation of a constraint which failed to converge, or the maximum violation
  ///        of any constraint within a group that failed to converge, for each failed RATTLE group
  ///        reported by getRattleFailures().
  ///
  /// \param tier  Obtain results from the CPU host or GPU device memory
  std::vector<float> getRattleViolations(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  
  /// \brief Get the violation of a constraint which failed to converge, or the maximum violation
  ///        of any constraint within a group that failed to converge, for each failed SHAKE group
  ///        reported by getRattleFailures().
  ///
  /// \param tier  Obtain results from the CPU host or GPU device memory
  std::vector<float> getShakeViolations(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the object's abstract.
  ///
  /// Overloaded:
  ///   - Obtain a read-only abstrat from a const object
  ///   - Obtain a mutable abstract from a non-const object
  ///
  /// \param tier  Indicate whether to target pointers to data on the CPU host or GPU device
  /// \{
  const WatcherReader data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  WatcherWriter data(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}
  
  /// \brief Set the threshold at which the object will report large forces.
  ///
  /// \param force_threshold_in  The magnitude of a force which will be reported as a significant
  ///                            event for subsequent error analysis
  void setForceThreshold(float force_threshold_in);

  /// \brief Set the threshold at which the object will report high particle speeds.
  ///
  /// \param speed_threshold_in  The magnitude of a particle velocity which will be reported as a
  ///                            significant event for subsequent error analysis
  void setSpeedThreshold(float speed_threshold_in);
  
private:

  // General parameters
  ExceptionResponse policy;  ///< The response to take in response to bad setup requests (may be
                             ///<   taken from the general program input policy)
  int max_reports;           ///< The maximum number of reports that the object will record
  
  // Track large forces on any given atom
  float force_threshold;           ///< The threshold for the magnitude of a force at which it will
                                   ///<   be reported as a significant event.
  Hybrid<int> large_force_count;   ///< The number of instances of large forces observed in the
                                   ///<   set of all simulations
  Hybrid<float4> large_forces;     ///< An array of large forces observed during all simulations,
                                   ///<   holding the Cartesian X, Y, and Z components of the force
                                   ///<   on some atom in the "x", "y", and "z" members of the
                                   ///<   tuple (as real-valued numbers).  The "w" member of the
                                   ///<   tuple holds the atom index cast bitwisae to float.
  Hybrid<int> large_force_steps;   ///< The step numbers at which each large force is observed
  Hybrid<int> large_force_stages;  ///< The points in the integration cycle at which each large
                                   ///<   force is encountered.

  // Track large velocities on any given atom
  float speed_threshold;             ///< The thresholld for the magnitude of a velocity at which
                                     ///<   it will be reported as a significant event
  Hybrid<int> high_speed_count;      ///< The number of instances of high velocities observed in
                                     ///<   the set of all simulations
  Hybrid<float4> high_speeds;        ///< An array of high velocities observed during all
                                     ///<   simulations, holding Cartesian X, Y, and Z components
                                     ///<   of the force on some atom in the "x", "y", and "z"
                                     ///<   members of teh tuple (as real-valued numbers).  The "w"
                                     ///<   member holds the atom index cast bitwise to float.
  Hybrid<int> high_speed_steps;      ///< The step numbers at which each high velocity is observed
  Hybrid<int> high_speed_stages;     ///< The points in the integration cycle at which each high
                                     ///<   velocity is encountered.
  
  // Track failed convergence in SHAKE or RATTLE groups
  Hybrid<int> failed_rattle_count;      ///< Total number of failed attempts to converge RATTLE
                                        ///<   groups
  Hybrid<int> failed_shake_count;       ///< Total number of failed attempts to converge SHAKE
                                        ///<   groups
  Hybrid<uint2> rattle_group_failures;  ///< The group indices (as given in the topology synthesis)
                                        ///<   of each constraint group that has failed to converge
                                        ///<   within the specified number of iterations in the "x"
                                        ///<   member of the tuple, followed by the step number
                                        ///<   upon which the attempt failed in the "y" member
  Hybrid<uint2> shake_group_failures;   ///< The group indices (as given in the topology synthesis)
                                        ///<   of each constraint group that has failed to converge
                                        ///<   within the specified number of iterations in the "x"
                                        ///<   member of the tuple, followed by the step number
                                        ///<   upon which the attempt failed in the "y" member
  Hybrid<float> rattle_violations;      ///< The maximum projection of any particle's velocity
                                        ///<   along the rigid bond axis for each of the RATTLE
                                        ///<   groups named in rattle_group_failures
  Hybrid<float> shake_violations;       ///< The maximum deviation from the target by which each of
                                        ///<   of the SHAKE groups named in shake_group_failures
                                        ///<   failed to converge

  // Track the amount of momentum purged with each cycle, if requested
  bool track_momentum_purge;
  Hybrid<float> x_velocity_purge;  ///< Translational velocity along the Cartesian X axis purged
                                   ///<   from each system, on each cycle
  Hybrid<float> y_velocity_purge;  ///< Translational velocity along the Cartesian Y axis purged
                                   ///<   from each system, on each cycle
  Hybrid<float> z_velocity_purge;  ///< Translational velocity along the Cartesian Z axis purged
                                   ///<   from each system, on each cycle
  Hybrid<float> x_angular_purge;   ///< Angular momentum about the X axis purged for each system
                                   ///<   on each cycle that the momentum purge is called.  This
                                   ///<   array is only allocated if the topology synthesis
                                   ///<   shows that the systems have no boundary conditions.
  Hybrid<float> y_angular_purge;   ///< Angular momentum about the Y axis purged for each system
                                   ///<   on each cycle that the momentum purge is called.
  Hybrid<float> z_angular_purge;   ///< Angular momentum about the Z axis purged for each system
                                   ///<   on each cycle that the momentum purge is called.

  /// ARRAY-kind Hybrid object targeted by various int-type Hybrids above.
  Hybrid<int> int_data;

  /// ARRAY-kind Hybrid object targeted by various float-type Hybrids above.
  Hybrid<float> float_data;

  /// Pointer to the associated coordinate synthesis
  PhaseSpaceSynthesis *poly_ps_ptr;

  /// Pointer to the associated topology synthesis
  AtomGraphSynthesis *poly_ag_ptr;
  
  /// \brief Allocate memory for the object.  Many of the arrays are POINTER-kind Hybrids
  ///        targeting the internal ARRAY-kind member variables.
  void allocate();

  /// \brief Validate the chosen force reporting threshold.
  void validateForceThreshold() const;

  /// \brief Validate the chosen speed reporting threshold.
  void validateSpeedThreshold() const;
};

} // namespace review
} // namespace stormm

#endif
