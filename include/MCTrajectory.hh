////////////////////////////////////////////////////////////////////////
/// \file  MCTrajectory.h
/// \version $Id: MCTrajectory.h,v 1.6 2012-11-01 19:18:11 brebel Exp $
/// \brief Trajectory class
///
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

/// This class describes the trajectory of a particle created in the
/// Monte Carlo simulation.  It generally behaves like a
/// vector< pair<TLorentzVector,TLorentzVector> >, where the first
/// TLorentzVector is the position and the seoond is the momentum,
/// with the following additions:

/// - Methods Position(int) and Momentum(int) for those who are unfamiliar with the
///   concept of "first" and "second" as used with STL pairs:
///      sim::Trajectory* trajectory = simb::MCParticle.Trajectory();
///      int numberOfPonts = trajectory->size();
///      for (int i=0; i<numberOfPoints; ++i)
///        {
///           TLorentzVector position = trajectory->Position(i);
///           TLorentzVector momentum = trajectory->Momentum(i);
///        }
///   The STL equivalent to the above statements (more efficient):
///      sim::Trajectory* trajectory = simb::MCParticle.Trajectory();
///      for ( sim::Trajectory::const_iterator i = trajectory->begin();
///            i != trajectory->end(); ++i )
///        {
///            const TLorentzVector& position = (*i).first;
///            const TLorentzVector& momentum = (*i).second;
///        }

/// - As above, but for each position or momentum component; e.g.,
///   trajectory->X(i).

/// - In addition to push_back(pair< TLorentzVector, TLorentzVector>),
///   there's also push_back(TLorentzVector,TLorentzVector) and
///   Add(TLorentzVector,TLorentzVector).  They all do the same thing:
///   add another point to the trajectory.

/// - Print() and operator<< methods for ROOT display and ease of
///   debugging.

/// There are no units defined in this class.  If it's used with
/// Geant4, the units will be (mm,ns,GeV), but this class does not
/// enforce this.

#ifndef SIMB_MCTRAJECTORY_H
#define SIMB_MCTRAJECTORY_H

#include <vector>
#include <iostream>

#include <TLorentzVector.h>

namespace simb {

  class MCTrajectory {
  public:
    /// Some type definitions to make life easier, and to help "hide"
    /// the implementation details.  (If you're not familiar with STL,
    /// you can ignore these definitions.)
    typedef std::vector< std::pair<TLorentzVector, TLorentzVector> >  list_type;
    typedef list_type::value_type                                     value_type;
    typedef list_type::iterator                                       iterator;
    typedef list_type::const_iterator                                 const_iterator;
    typedef list_type::reverse_iterator                               reverse_iterator;
    typedef list_type::const_reverse_iterator                         const_reverse_iterator;
    typedef list_type::size_type                                      size_type;
    typedef list_type::difference_type                                difference_type;
    typedef std::vector< std::pair<size_t, unsigned char> >           ProcessMap;
    /// Standard constructor: Start with initial position and momentum
    /// of the particle.
    MCTrajectory();

  private:
    list_type  ftrajectory;        ///< The list of trajectory points
    ProcessMap fTrajectoryProcess; ///< map of the scattering process to index
                                   ///< in ftrajectory for a given point

  public:

    MCTrajectory( const TLorentzVector& vertex,
                  const TLorentzVector& momentum );

    /// The accessor methods described above.
    const TLorentzVector& Position( const size_type ) const;
    const TLorentzVector& Momentum( const size_type ) const;
    double  X( const size_type i ) const;
    double  Y( const size_type i ) const;
    double  Z( const size_type i ) const;
    double  T( const size_type i ) const;
    double Px( const size_type i ) const;
    double Py( const size_type i ) const;
    double Pz( const size_type i ) const;
    double  E( const size_type i ) const;

    double TotalLength() const;

    friend std::ostream& operator<< ( std::ostream& output, const MCTrajectory& );

    /// Standard STL methods, to make this class look like an STL map.
    /// Again, if you don't know STL, you can just ignore these
    /// methods.
    iterator               begin();
    const_iterator         begin()      const;
    iterator               end();
    const_iterator         end()        const;
    reverse_iterator       rbegin();
    const_reverse_iterator rbegin()     const;
    reverse_iterator       rend();
    const_reverse_iterator rend()       const;

    size_type size()                    const;
    bool      empty()                   const;
    void      swap(simb::MCTrajectory& other);
    void      clear();

    // Note that there's no non-const version of operator[] or at() here; once
    // you've added a point to a trajectory, you can't modify it.
    const value_type& operator[](const size_type i) const;
    const value_type& at(const size_type i)         const;

    /// The only "set" methods for this class; once you've added a
    /// trajectory point, you can't take it back.
    void push_back(value_type const& v );
    void push_back(TLorentzVector const& p,
                   TLorentzVector const& m );
    void Add(TLorentzVector const& p,
             TLorentzVector const& m );
    void Add(TLorentzVector const& p,
             TLorentzVector const& m,
             std::string    const& process,
             bool keepTransportation = false);

    unsigned char        ProcessToKey(std::string   const& process) const;
    std::string          KeyToProcess(unsigned char const& key)     const;
    ProcessMap    const& TrajectoryProcesses()                      const;

    /// Remove points from trajectory. Straight line interpolation between the
    /// remaining points will pass no further than \a margin from removed
    /// points.
    void Sparsify(double margin = .1, bool keep_second_to_last = false);

  };

} // namespace simb

inline double           simb::MCTrajectory::X ( const size_type i ) const { return Position(i).X();      }
inline double           simb::MCTrajectory::Y ( const size_type i ) const { return Position(i).Y();      }
inline double           simb::MCTrajectory::Z ( const size_type i ) const { return Position(i).Z();      }
inline double           simb::MCTrajectory::T ( const size_type i ) const { return Position(i).T();      }
inline double           simb::MCTrajectory::Px( const size_type i ) const { return Momentum(i).Px();     }
inline double           simb::MCTrajectory::Py( const size_type i ) const { return Momentum(i).Py();     }
inline double           simb::MCTrajectory::Pz( const size_type i ) const { return Momentum(i).Pz();     }
inline double           simb::MCTrajectory::E ( const size_type i ) const { return Momentum(i).E();      }

inline simb::MCTrajectory::iterator               simb::MCTrajectory::begin()             { return ftrajectory.begin();  }
inline simb::MCTrajectory::const_iterator         simb::MCTrajectory::begin()       const { return ftrajectory.begin();  }
inline simb::MCTrajectory::iterator               simb::MCTrajectory::end()               { return ftrajectory.end();    }
inline simb::MCTrajectory::const_iterator         simb::MCTrajectory::end()         const { return ftrajectory.end();    }
inline simb::MCTrajectory::reverse_iterator       simb::MCTrajectory::rbegin()            { return ftrajectory.rbegin(); }
inline simb::MCTrajectory::const_reverse_iterator simb::MCTrajectory::rbegin()      const { return ftrajectory.rbegin(); }
inline simb::MCTrajectory::reverse_iterator       simb::MCTrajectory::rend()              { return ftrajectory.rend();   }
inline simb::MCTrajectory::const_reverse_iterator simb::MCTrajectory::rend()        const { return ftrajectory.rend();   }
inline simb::MCTrajectory::size_type              simb::MCTrajectory::size()        const { return ftrajectory.size();   }
inline bool                                       simb::MCTrajectory::empty()       const { return ftrajectory.empty();  }
inline void                                       simb::MCTrajectory::clear()             { ftrajectory.clear();         }
inline void                                       simb::MCTrajectory::swap(simb::MCTrajectory& other)
                                                                                { ftrajectory.swap( other.ftrajectory ); }

inline const simb::MCTrajectory::value_type&      simb::MCTrajectory::operator[](const simb::MCTrajectory::size_type i) const
                                                                                                { return ftrajectory[i]; }

inline const simb::MCTrajectory::value_type&      simb::MCTrajectory::at(const simb::MCTrajectory::size_type i)         const
                                                                                             { return ftrajectory.at(i); }

inline void                                       simb::MCTrajectory::push_back(const simb::MCTrajectory::value_type& v )
                                                                                             { ftrajectory.push_back(v); }

inline void                                       simb::MCTrajectory::push_back(const TLorentzVector& p,
                                                                                const TLorentzVector& m )
                                                         { ftrajectory.push_back( simb::MCTrajectory::value_type(p,m) ); }

inline void                                       simb::MCTrajectory::Add(const TLorentzVector& p,
                                                                          const TLorentzVector& m )    { push_back(p,m); }

inline simb::MCTrajectory::ProcessMap    const&   simb::MCTrajectory::TrajectoryProcesses() const { return fTrajectoryProcess; }

#endif // SIMB_MCTRAJECTORY_H
