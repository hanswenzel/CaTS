////////////////////////////////////////////////////////////////////////
/// \file  MCTrajectory.cxx
/// \brief Container of trajectory information for a particle
///
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

#include "MCTrajectory.hh"

#include <TLorentzVector.h>

#include <stdexcept>
#include <sstream>
#include <cmath>
#include <deque>
#include <iterator>
#include <vector>
#include <set>
#include <map>

namespace simb
{

  // Nothing special need be done for the default constructor or destructor.
  MCTrajectory::MCTrajectory()
    : ftrajectory()
  {}

  //----------------------------------------------------------------------------
  MCTrajectory::MCTrajectory(const TLorentzVector& position, const TLorentzVector& momentum)
  {
    ftrajectory.push_back(value_type(position, momentum));
  }

  //----------------------------------------------------------------------------
  const TLorentzVector& MCTrajectory::Position(const size_type index) const
  {
    const_iterator i = ftrajectory.begin();
    std::advance(i, index);
    return (*i).first;
  }

  //----------------------------------------------------------------------------
  const TLorentzVector& MCTrajectory::Momentum(const size_type index) const
  {
    const_iterator i = ftrajectory.begin();
    std::advance(i, index);
    return (*i).second;
  }

  //----------------------------------------------------------------------------
  double MCTrajectory::TotalLength() const
  {
    const int N = size();
    if(N < 2)
      return 0;

    // We take the sum of the straight lines between the trajectory points
    double dist = 0;
    for(int n = 0; n < N - 1; ++n)
    {
      dist += (Position(n + 1) - Position(n)).Vect().Mag();
    }

    return dist;
  }

  //----------------------------------------------------------------------------
  std::ostream& operator<<(std::ostream& output, const MCTrajectory& list)
  {
    // Determine a field width for the voxel number.
    MCTrajectory::size_type numberOfTrajectories = list.size();
    int numberOfDigits = (int) std::log10((double) numberOfTrajectories) + 1;

    // A simple header.
    output.width(numberOfDigits);
    output << "#"
           << ": < position (x,y,z,t), momentum (Px,Py,Pz,E) >" << std::endl;

    // Write each trajectory point on a separate line.
    MCTrajectory::size_type nTrajectory = 0;
    for(MCTrajectory::const_iterator trajectory = list.begin(); trajectory != list.end();
        ++trajectory, ++nTrajectory)
    {
      output.width(numberOfDigits);
      output << nTrajectory << ": "
             << "< (" << (*trajectory).first.X() << "," << (*trajectory).first.Y() << ","
             << (*trajectory).first.Z() << "," << (*trajectory).first.T() << ") , ("
             << (*trajectory).second.Px() << "," << (*trajectory).second.Py() << ","
             << (*trajectory).second.Pz() << "," << (*trajectory).second.E() << ") >" << std::endl;
    }

    return output;
  }

  //----------------------------------------------------------------------------
  unsigned char MCTrajectory::ProcessToKey(std::string const& process) const
  {
    unsigned char key = 0;

    if(process.compare("hadElastic") == 0)
      key = 1;
    else if(process.compare("pi-Inelastic") == 0)
      key = 2;
    else if(process.compare("pi+Inelastic") == 0)
      key = 3;
    else if(process.compare("kaon-Inelastic") == 0)
      key = 4;
    else if(process.compare("kaon+Inelastic") == 0)
      key = 5;
    else if(process.compare("protonInelastic") == 0)
      key = 6;
    else if(process.compare("neutronInelastic") == 0)
      key = 7;
    else if(process.compare("CoulombScat") == 0)
      key = 8;
    else if(process.compare("nCapture") == 0)
      key = 9;
    else if(process.compare("Transportation") == 0)
      key = 10;

    return key;
  }

  //----------------------------------------------------------------------------
  std::string MCTrajectory::KeyToProcess(unsigned char const& key) const
  {
    std::string process("Unknown");

    if(key == 1)
      process = "hadElastic";
    else if(key == 2)
      process = "pi-Inelastic";
    else if(key == 3)
      process = "pi+Inelastic";
    else if(key == 4)
      process = "kaon-Inelastic";
    else if(key == 5)
      process = "kaon+Inelastic";
    else if(key == 6)
      process = "protonInelastic";
    else if(key == 7)
      process = "neutronInelastic";
    else if(key == 8)
      process = "CoulombScat";
    else if(key == 9)
      process = "nCapture";
    else if(key == 10)
      process = "Transportation";

    return process;
  }

  //----------------------------------------------------------------------------
  void MCTrajectory::Add(TLorentzVector const& p, TLorentzVector const& m,
                         std::string const& process, bool keepTransportation)
  {
    // add the the momentum and position, then get the location of the added
    // bits to store the process
    this->push_back(p, m);

    size_t insertLoc = ftrajectory.size() - 1;

    auto key = this->ProcessToKey(process);

    // only add a process to the list if it is defined, ie one of the values
    // allowed in the ProcessToKey() method
    //
    // Also, keep 10 (transportation) if the flag allows
    if(key > 0 && (key != 10 || keepTransportation))
      fTrajectoryProcess.push_back(std::make_pair(insertLoc, key));

    return;
  }

  //----------------------------------------------------------------------------
  void MCTrajectory::Sparsify(double margin, bool keep_second_to_last)
  {
    // This is a divide-and-conquer algorithm. If the straight line between two
    // points is close enough to all the intermediate points, then just keep
    // the endpoints. Otherwise, divide the range in two and try again.

    // We keep the ranges that need checking in "toCheck". If a range is good
    // as-is, we put just the starting point in "done". The end-point will be
    // taken care of by the next range.

    // Need at least three points to think of removing one
    // D.R. -- let's up this to four points before we start removing
    //      -- this is helpful when retrieving the energy of the particle prior
    //      -- to a final interaction : (Start, p1, ..., p_(n-1), End)
    if(size() <= 3 && keep_second_to_last)
      return;
    else if(size() <= 2)
      return;

    // Deal in terms of distance-squared to save some sqrts
    margin *= margin;

    // Deque because we add things still to check on the end, and pop things
    // we've checked from the front.
    std::deque<std::pair<int, int>> toCheck;
    // Start off by trying to replace the whole trajectory with just the
    // endpoints.
    toCheck.push_back(std::make_pair(0, size() - 1));

    // use a std::set to keep track of which indices of points we want to
    // keep because the set does not allow duplicates and it keeps items in
    // order
    std::set<int> done;
    if(keep_second_to_last)
      done.insert(size() - 2);  // -- D.R. store the penultimate point

    while(!toCheck.empty())
    {
      const int loIdx = toCheck.front().first;
      const int hiIdx = toCheck.front().second;
      toCheck.pop_front();

      // Should never have been given a degenerate range
      if(hiIdx < loIdx + 2)
        throw std::runtime_error("MCTrajectory: Degnerate range in Sparsify method");
      const TVector3 loVec = at(loIdx).first.Vect();
      const TVector3 hiVec = at(hiIdx).first.Vect();

      const TVector3 dir = (hiVec - loVec).Unit();

      // Are all the points in between close enough?
      bool ok = true;
      for(int i = loIdx + 1; i < hiIdx; ++i)
      {
        const TVector3 toHere = at(i).first.Vect() - loVec;
        // Perpendicular distance^2 from the line joining loVec to hiVec
        const double impact = (toHere - dir.Dot(toHere) * dir).Mag2();
        if(impact > margin)
        {
          ok = false;
          break;
        }
      }

      if(ok)
      {
        // These points adequately represent this range
        done.insert(loIdx);
      }
      else
      {
        // Split in half
        const int midIdx = (loIdx + hiIdx) / 2;
        // Should never have a range this small
        if(midIdx == loIdx)
          throw std::runtime_error("MCTrajectory: Midpoint in sparsification is same as lowpoint");
        if(midIdx == hiIdx)
          throw std::runtime_error("MCTrajectory: Midpoint in sparsification is same as hipoint");

        // The range can be small enough that upon splitting, the new ranges
        // are degenerate, and should just be written out straight away. Check
        // for those cases.

        if(midIdx == loIdx + 1)
        {
          done.insert(loIdx);
        }
        else
        {
          toCheck.push_back(std::make_pair(loIdx, midIdx));
        }

        if(midIdx == hiIdx - 1)
        {
          done.insert(midIdx);
        }
        else
        {
          toCheck.push_back(std::make_pair(midIdx, hiIdx));
        }
      }
    }  // end while

    // now make sure we have not left out any of the indices where interesting
    // processes were recorded
    std::map<size_t, unsigned char> processMap;
    for(auto itr : fTrajectoryProcess)
    {
      done.insert(itr.first);
      processMap[itr.first] = itr.second;
    }

    // Look up the trajectory points at the stored indices, write them to a new
    // trajectory
    const unsigned int I = done.size();
    list_type newTraj;
    newTraj.reserve(I + 1);

    // make a new process map as well based on these points
    ProcessMap newProcMap;

    for(auto idx : done)
    {
      newTraj.push_back(at(idx));
      if(processMap.count(idx) > 0)
      {
        newProcMap.push_back(std::make_pair(newTraj.size() - 1, processMap.find(idx)->second));
      }
    }

    // Remember to add the very last point in if it hasn't already been added
    if(done.count(ftrajectory.size() - 1) < 1)
      newTraj.push_back(*rbegin());

    // Replace trajectory and fTrajectoryProcess with new versions
    std::swap(ftrajectory, newTraj);
    std::swap(fTrajectoryProcess, newProcMap);

    return;
  }

}  // namespace simb
