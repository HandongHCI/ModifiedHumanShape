/**
    This file is part of the implementation of the 3D human shape model as described in the paper:

    Leonid Pishchulin, Stefanie Wuhrer, Thomas Helten, Christian Theobalt and Bernt Schiele
    Building Statistical Shape Spaces for 3D Human Modeling
    ArXiv, March 2015

    Please cite the paper if you are using this code in your work.
    
    Contributors: Arjun Jain, Juergen Gall, Leonid Pishchulin.
    This code also uses open source functionality for matrix operations by Thomas Brox.

    The code may be used free of charge for non-commercial and
    educational purposes, the only requirement is that this text is
    preserved within the derivative work. For any other purpose you
    must contact the authors for permission. This code may not be
    redistributed without permission from the authors.
*/

#include "NMath.h"
#include "NRBM.h"
#include "CTMesh.h"
#include <assert.h>
#include <set>
#include <map>
#include <string.h>

std::vector<CMatrix<double> > CMesh::eigenVectors;
CMatrix<float> CMesh::weightMatrix;

using namespace std;
// C J O I N T -----------------------------------------------------------------

// constructor
CJoint::CJoint(CVector<float>& aDirection, CVector<float>& aPoint, int aParent) {
	mDirection = aDirection;
	mPoint = aPoint;
	mMoment = aPoint/aDirection;
	mParent = aParent;
}

// rigidMotion
void CJoint::rigidMotion(CMatrix<float>& M) {
	CMatrix<float> R(3,3);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			R(i,j) = M(i,j);
	mDirection = R*mDirection;
	mPoint = R*mPoint;
	mPoint(0) += M(3,0); mPoint(1) += M(3,1); mPoint(2) += M(3,2);
	mMoment = mPoint/mDirection;
}

// angleToMatrix
void CJoint::angleToMatrix(float aAngle, CMatrix<float>& M) {
	CMatrix<float> omegaHat(3,3);
	omegaHat.data()[0] = 0.0;            omegaHat.data()[1] = -mDirection(2); omegaHat.data()[2] = mDirection(1);
	omegaHat.data()[3] = mDirection(2);  omegaHat.data()[4] = 0.0;            omegaHat.data()[5] = -mDirection(0);
	omegaHat.data()[6] = -mDirection(1); omegaHat.data()[7] = mDirection(0);  omegaHat.data()[8] = 0.0;
	CMatrix<float> omegaT(3,3);
	for (int j = 0; j < 3; j++)
	  for (int i = 0; i < 3; i++)
	    omegaT(i,j) = mDirection(i)*mDirection(j);
	CMatrix<float> R(3,3);
	R = (omegaHat*(float)sin(aAngle))+((omegaHat*omegaHat)*(float)(1.0-cos(aAngle)));
	R(0,0) += 1.0; R(1,1) += 1.0; R(2,2) += 1.0;
	CMatrix<float> T(R); T *= -1.0;
	T(0,0) += 1.0; T(1,1) += 1.0; T(2,2) += 1.0;
	CVector<float> t(3);
	t = T*(mDirection/mMoment)+((omegaT*mMoment)*aAngle);
	M.data()[0] = R.data()[0]; M.data()[1] = R.data()[1]; M.data()[2] = R.data()[2]; M.data()[3] = t(0);
	M.data()[4] = R.data()[3]; M.data()[5] = R.data()[4]; M.data()[6] = R.data()[5]; M.data()[7] = t(1);
	M.data()[8] = R.data()[6]; M.data()[9] = R.data()[7]; M.data()[10]= R.data()[8]; M.data()[11]= t(2);
	M.data()[12]= 0.0;         M.data()[13]= 0.0;         M.data()[14]= 0.0;         M.data()[15]= 1.0;
}

// operator =
CJoint& CJoint::operator=(CJoint& aCopyFrom) {
	mDirection = aCopyFrom.mDirection;
	mPoint = aCopyFrom.mPoint;
	mMoment = aCopyFrom.mMoment;
	mParent = aCopyFrom.mParent;
	return *this;
}

// C M E S H M O T I O N -------------------------------------------------------

// reset
void CMeshMotion::reset(int aJointNumber) {
	mRBM.setSize(4,4);
	mRBM = 0.0;
	mRBM(0,0) = 1.0; mRBM(1,1) = 1.0; mRBM(2,2) = 1.0; mRBM(3,3) = 1.0;
	mPoseParameters.setSize(6+aJointNumber);
	mPoseParameters = 0.0;
}

// print
void CMeshMotion::print() {
  std::cout << mRBM << std::endl;
  std::cout << mPoseParameters;
}

// writeToFile
void CMeshMotion::writeToFile() {
	//char buffer[200];
	//sprintf(buffer,"%s%03d.txt",(NShow::mResultDir+"RBM").c_str(),NShow::mImageNo);
	//mRBM.writeToTXT(buffer);
	//sprintf(buffer,"%s%03d.txt",(NShow::mResultDir+"Params").c_str(),NShow::mImageNo);
	//mPoseParameters.writeToTXT(buffer);
}

// operator =
CMeshMotion& CMeshMotion::operator=(const CMeshMotion& aCopyFrom) {
	mRBM = aCopyFrom.mRBM;
	mPoseParameters = aCopyFrom.mPoseParameters;
	return *this;
}

// C M E S H -------------------------------------------------------------------

// constructor
CMesh::CMesh(int aPatches, int aPoints) {
	mPoints.resize(aPoints);
	mPatch.resize(aPatches);
	mNumPoints=aPoints;
	mNumPatch=aPatches;

	for(int i=0;i<aPoints;i++)
	{
		mPoints[i].setSize(4);
	}
	for(int i=0;i<aPatches;i++)
	{
		mPatch[i].setSize(3);
	}
}

// destructor
CMesh::~CMesh() { 
	//cout << "DEL-CMESH...";
	//if(mPointsS !=0) delete[] mPointsS;
	//if(mPoints !=0) delete[] mPoints;
	//if(mPatch !=0) delete[] mPatch;
	//if(mAppearanceField !=0) delete[] mAppearanceField;
	//cout << "...ok\n"; 
}

// readModel
bool CMesh::readModelHybOpt(const char* aFilename, bool smooth) {
#if 0
#ifndef EXEC_FAST
	cout << "Read Model... ";
#endif

	std::ifstream aStream(aFilename);
	if(!aStream.is_open()) {
		cerr << "Could not open " << aFilename << endl;
		return false;
	}

	int aXSize, aYSize, size = 4;

	CJoint Joint;

	aStream >> mNumPoints; 
	aStream >> mNumPatch; 
	aStream >> mJointNumber;
	if(smooth) {
		aStream >> mNumSmooth;
		if(mNumSmooth == 1)
			smooth = false;
		else
			size = 3 + mNumSmooth * 2;
	} else
		mNumSmooth = 1;

	weightMatrix.setSize(mNumPoints, mJointNumber);
	weightMatrix = 0;

#ifndef EXEC_FAST
	cout << mNumPoints << " " << mNumPatch << " " << mJointNumber << " " << mNumSmooth << endl;
#endif

	CVector<bool> BoundJoints(mJointNumber+1,false);

	mJoint.setSize(mJointNumber+1);

	// Read mesh components
	mPoints .resize(mNumPoints); 
	mPatch.resize(mNumPatch);

	for(int i=0;i<mNumPoints;i++) {
		mPoints[i].setSize(size);
	} 

	for(int i=0;i<mNumPatch;i++) {
		mPatch[i].setSize(3);
	}

	for (int i = 0; i < mNumPoints; i++) {
		aStream >> mPoints[i][0]; 
		aStream >> mPoints[i][1];
		aStream >> mPoints[i][2]; 
		if(smooth) {
			for(int n = 0; n < mNumSmooth * 2; n++)
				aStream >> mPoints[i][3 + n];
		} else
			aStream >> mPoints[i][3]; 
		BoundJoints((int)mPoints[i][3]) = true;
	}

	for(int i0 = 0; i0 < mNumPoints; i0++)
	{
		for(int i1 = 3; i1 < mNumSmooth * 2 + 3; i1 += 2)
		{
			int jointID = int(mPoints[i0][i1]);
			float weight = mPoints[i0][i1 + 1];
			weightMatrix(i0, jointID - 1) += weight;
		}
	} 

	for (int i = 0; i < mNumPatch; i++) {
		aStream >> mPatch[i][0];
		aStream >> mPatch[i][1];
		aStream >> mPatch[i][2];
	}

	// set bounds
	int count = 0;
	mJointMap.setSize(mJointNumber+1);
	mJointMap = -1;
	for(int j = 0; j<=mJointNumber; ++j)
		if(BoundJoints(j))
			mJointMap(j) = count++;     

	mBounds.setSize(count,9,3);
	mNeighbor.setSize(count, count); 
	CMatrix<float> minV(count,3,100000);
	CMatrix<float> maxV(count,3,-100000);

	for (int i = 0; i < mNumPoints; i++) {
		int index = mJointMap((int)mPoints[i][3]);

		if ( mPoints[i][0]<minV(index,0) ) minV(index,0) = mPoints[i][0];
		if ( mPoints[i][1]<minV(index,1) ) minV(index,1) = mPoints[i][1];
		if ( mPoints[i][2]<minV(index,2) ) minV(index,2) = mPoints[i][2];
		if ( mPoints[i][0]>maxV(index,0) ) maxV(index,0) = mPoints[i][0];
		if ( mPoints[i][1]>maxV(index,1) ) maxV(index,1) = mPoints[i][1];
		if ( mPoints[i][2]>maxV(index,2) ) maxV(index,2) = mPoints[i][2];
	}

	for(int i=0; i<count; ++i) {
		mBounds(i,0,0) = mBounds(i,1,0) = mBounds(i,2,0) = mBounds(i,3,0) = minV(i,0);
		mBounds(i,4,0) = mBounds(i,5,0) = mBounds(i,6,0) = mBounds(i,7,0) = maxV(i,0);
		mBounds(i,0,1) = mBounds(i,1,1) = mBounds(i,4,1) = mBounds(i,5,1) = minV(i,1);
		mBounds(i,2,1) = mBounds(i,3,1) = mBounds(i,6,1) = mBounds(i,7,1) = maxV(i,1);
		mBounds(i,0,2) = mBounds(i,2,2) = mBounds(i,4,2) = mBounds(i,6,2) = minV(i,2);
		mBounds(i,1,2) = mBounds(i,3,2) = mBounds(i,5,2) = mBounds(i,7,2) = maxV(i,2);
	}

	mCenter.setSize(4);
	mCenter[0] = (maxV(0,0) + minV(0,0))/2.0f;
	mCenter[1] = (maxV(0,1) + minV(0,1))/2.0f;
	mCenter[2] = (maxV(0,2) + minV(0,2))/2.0f;
	mCenter[3] = 1.0f;

	for(int j=0; j<mJointMap.size();++j)
		if(mJointMap(j)>=0)
			mBounds(mJointMap(j),8,0) = j;

	// Read joints
	int dummy;
	CVector<float> aDirection(3);
	CVector<float> aPoint(3);
	for (int aJointID = 1; aJointID <= mJointNumber; aJointID++) {
		aStream >> dummy; // ID
		aStream >> aDirection(0) >> aDirection(1) >> aDirection(2);
		aStream >> aPoint(0) >> aPoint(1) >> aPoint(2);
		mJoint(aJointID).set(aDirection,aPoint);
		aStream >> mJoint(aJointID).mParent;
	}
	// Determine which joint motion is influenced by parent joints
	mInfluencedBy.setSize(mJointNumber+1,mJointNumber+1);
	mInfluencedBy = false;
	for (int j = 0; j <= mJointNumber; j++)
		for (int i = 0; i <= mJointNumber; i++) {
			if (i == j) mInfluencedBy(i,j) = true;
			if (isParentOf(j,i)) mInfluencedBy(i,j) = true;
		}

		mNeighbor = false;
		for (int i = 0; i < mNeighbor.xSize(); i++) {
			mNeighbor(i,i) = true;
			int jID = (int)mBounds(i,8,0);
			for (int j = jID-1; j >= 0; --j) {
				if(mJointMap(j)>=0 && mInfluencedBy(jID,j)) {
					mNeighbor(i,mJointMap(j)) = true;
					mNeighbor(mJointMap(j),i) = true;
					break;
				}
			}
		}

		mEndJoint.setSize(mJointNumber+1);
		mEndJoint.fill(true);
		for(int i=1; i<=mJointNumber; ++i)
			mEndJoint[mJoint(i).mParent]=false;

		for(int i=1; i<=mJointNumber; ++i)
			if(mEndJoint[i]==true) {
				int j = i;
				while(mJointMap[--j]==-1) {
					mEndJoint[j]=true;
				}
			}

			cout << "End Joints:" << endl;
			cout << mEndJoint << endl;

			mAccumulatedMotion.reset(mJointNumber);
			mCurrentMotion.reset(mJointNumber);

#ifndef EXEC_FAST
			cout << mNeighbor << endl;
			cout << "ok" << endl;
#endif
#endif
			return true;
}

// readModel
bool CMesh::readModel(const char* aFilename, bool smooth) 
{
#ifndef EXEC_FAST
	cout << "Read Model... ";
#endif

	std::ifstream aStream(aFilename);
	if(!aStream.is_open()) {
		cerr << "Could not open " << aFilename << endl;
		return false;
	}

	//int aXSize, aYSize
	int size = 4;

	CJoint Joint;

	aStream >> mNumPoints; 
	aStream >> mNumPatch; 
	aStream >> mJointNumber;
	if(smooth) {
		aStream >> mNumSmooth;
		if(mNumSmooth == 1)
			smooth = false;
		else
			size = 3 + mNumSmooth * 2;
	} else
		mNumSmooth = 1;

	weightMatrix.setSize(mNumPoints, mJointNumber);
	weightMatrix = 0;

#ifndef EXEC_FAST
	cout << mNumPoints << " " << mNumPatch << " " << mJointNumber << " " << mNumSmooth << endl;
#endif

	CVector<bool> BoundJoints(mJointNumber+1,false);
	mBoundJoints.setSize(mJointNumber+1);

	//mHistory.setSize(5);
	mJoint.setSize(mJointNumber+1);
	// Read mesh components
	mPoints.resize(mNumPoints); 
	mPatch.resize(mNumPatch);

	for(int i=0;i<mNumPoints;i++) {
		mPoints[i].setSize(size);
	} 

	for(int i=0;i<mNumPatch;i++) {
		mPatch[i].setSize(3);
	}

	//   mPointsS=new CVector<float>[mNumPoints];
	//   for(int i=0;i<mNumPoints;i++)
	//     {
	//       mPointsS[i].setSize(4);
	//     }
	for (int i = 0; i < mNumPoints; i++) {
		aStream >> mPoints[i][0]; 
		aStream >> mPoints[i][1];
		aStream >> mPoints[i][2]; 
		//~aj: FIX FIX
		for(int m0 = 0; m0 < 3; m0++)
			mPoints[i][m0] *= SKEL_SCALE_FACTOR;
		
		if(smooth) {
			for(int n = 0; n < mNumSmooth * 2; n++)
				aStream >> mPoints[i][3 + n];
		} else
			aStream >> mPoints[i][3]; 
		BoundJoints((int)mPoints[i][3]) = true;
	}
	for(int i0 = 0; i0 < mNumPoints; i0++)
	{
		for(int i1 = 3; i1 < mNumSmooth * 2 + 3; i1 += 2)
		{
			int jointID = int(mPoints[i0][i1]);
			float weight = mPoints[i0][i1 + 1];
			weightMatrix(i0, jointID - 1) += weight;
		}
	} 
	for (int i = 0; i < mNumPatch; i++) {
		aStream >> mPatch[i][0];
		aStream >> mPatch[i][1];
		aStream >> mPatch[i][2];
	}

	// set bounds
	int count = 0;
	mJointMap.setSize(mJointNumber+1);
	mJointMap = -1;
	for(int j = 0; j<=mJointNumber; ++j)
		if(BoundJoints(j))
			mJointMap(j) = count++;      
	
	mNoOfBodyParts = count;

#ifndef EXEC_FAST
	cout << "Bodyparts: " << count << endl;
	cout<<"\nmJointMap: "<<mJointMap;
	cout<<"\nBoundJoints: "<<BoundJoints;
#endif	
	mBoundJoints = BoundJoints;


	mBounds.setSize(count,9,3);
	mNeighbor.setSize(count, count); 
	CMatrix<float> minV(count,3,100000);
	CMatrix<float> maxV(count,3,-100000);

	mCenter.setSize(4);
	mCenter[0] = 0;
	mCenter[1] = 0;
	mCenter[2] = 0;
	mCenter[3] = 1.0f;

	for (int i = 0; i < mNumPoints; i++) {
		int index = mJointMap((int)mPoints[i][3]);

		if ( mPoints[i][0]<minV(index,0) ) minV(index,0) = mPoints[i][0];
		if ( mPoints[i][1]<minV(index,1) ) minV(index,1) = mPoints[i][1];
		if ( mPoints[i][2]<minV(index,2) ) minV(index,2) = mPoints[i][2];
		if ( mPoints[i][0]>maxV(index,0) ) maxV(index,0) = mPoints[i][0];
		if ( mPoints[i][1]>maxV(index,1) ) maxV(index,1) = mPoints[i][1];
		if ( mPoints[i][2]>maxV(index,2) ) maxV(index,2) = mPoints[i][2];

		mCenter[0]+=mPoints[i][0];
		mCenter[1]+=mPoints[i][1];
		mCenter[2]+=mPoints[i][2];

	}

	mCenter[0] /= (float)mNumPoints;
	mCenter[1] /= (float)mNumPoints;
	mCenter[2] /= (float)mNumPoints;
	
#ifndef EXEC_FAST
	cout << "mCenter="  << mCenter << endl;
#endif
	for(int i=0; i<count; ++i) {
		mBounds(i,0,0) = mBounds(i,1,0) = mBounds(i,2,0) = mBounds(i,3,0) = minV(i,0);
		mBounds(i,4,0) = mBounds(i,5,0) = mBounds(i,6,0) = mBounds(i,7,0) = maxV(i,0);
		mBounds(i,0,1) = mBounds(i,1,1) = mBounds(i,4,1) = mBounds(i,5,1) = minV(i,1);
		mBounds(i,2,1) = mBounds(i,3,1) = mBounds(i,6,1) = mBounds(i,7,1) = maxV(i,1);
		mBounds(i,0,2) = mBounds(i,2,2) = mBounds(i,4,2) = mBounds(i,6,2) = minV(i,2);
		mBounds(i,1,2) = mBounds(i,3,2) = mBounds(i,5,2) = mBounds(i,7,2) = maxV(i,2);
	}

	for(int j=0; j<mJointMap.size();++j)
		if(mJointMap(j)>=0)
			mBounds(mJointMap(j),8,0) = (float)j;

	// Read joints
	int dummy;
	CVector<float> aDirection(3);
	CVector<float> aPoint(3);
	
	for (int aJointID = 1; aJointID <= mJointNumber; aJointID++) {
		aStream >> dummy; // ID
		aStream >> aDirection(0) >> aDirection(1) >> aDirection(2);
		aStream >> aPoint(0) >> aPoint(1) >> aPoint(2);

		for(int m0 = 0; m0 < 3; m0++){
		  aPoint(m0) *= SKEL_SCALE_FACTOR;
			
		}
		mJoint(aJointID).set(aDirection,aPoint);
		aStream >> mJoint(aJointID).mParent;
	}
	// Determine which joint motion is influenced by parent joints
	mInfluencedBy.setSize(mJointNumber+1,mJointNumber+1);
	mInfluencedBy = false;
	for (int j = 0; j <= mJointNumber; j++)
	  for (int i = 0; i <= mJointNumber; i++) {
	    if (i == j) mInfluencedBy(i,j) = true;
	    if (isParentOf(j,i)) mInfluencedBy(i,j) = true;
	  }

	mNeighbor = false;
		for (int i = 0; i < mNeighbor.xSize(); i++) {
			mNeighbor(i,i) = true;
			int jID = (int)mBounds(i,8,0);
			for (int j = jID-1; j >= 0; --j) {
				if(mJointMap(j)>=0 && mInfluencedBy(jID,j)) {
					mNeighbor(i,mJointMap(j)) = true;
					mNeighbor(mJointMap(j),i) = true;
					break;
				}
			}
		}

		mEndJoint.setSize(mJointNumber+1);
		mEndJoint.fill(true);
		for(int i=1; i<=mJointNumber; ++i)
			mEndJoint[mJoint(i).mParent]=false;

		for(int i=1; i<=mJointNumber; ++i)
			if(mEndJoint[i]==true) {
				int j = i;
				while(mJointMap[--j]==-1 && j > 0) {
					mEndJoint[j]=true;
				}
			}
#ifndef EXEC_FAST
		cout << "End Joints:" << endl;
		cout << mEndJoint << endl;
#endif
		mAccumulatedMotion.reset(mJointNumber);
		mCurrentMotion.reset(mJointNumber);
		
		mCovered.setSize(mBounds.xSize());
		mExtremity.setSize(mBounds.xSize());
		for(int i=0; i<mExtremity.size(); ++i)
		  aStream >> mExtremity[i];
		for(int i=0; i<mCovered.size(); ++i)
		  aStream >> mCovered[i];

#ifndef EXEC_FAST
		cout << mNeighbor << endl;
		cout << "ok" << endl;
#endif


			return true;
}

void CMesh::centerModel() {

#ifndef EXEC_FAST
	cout << "RESET CENTER!\n";
#endif

	CVector<float> trans(4);
	trans(0) = -mCenter(0); trans(1) = -mCenter(1); trans(2) = -mCenter(2); trans(3) = 0;

	#ifndef EXEC_FAST
	cout << "trans:" << trans;
#endif

	for (int i = 0; i < mNumPoints; i++) {
		mPoints[i](0) += trans(0);
		mPoints[i](1) += trans(1);
		mPoints[i](2) += trans(2);
	}

	for (int i = 1; i <= mJointNumber; i++) {
		CVector<float> jPos(mJoint[i].getPoint());
		for(int j = 0; j < 3; ++j) 
		  jPos(j) += trans(j);
		mJoint[i].setPoint(jPos);
	}
	
	for(int i = 0; i < mBounds.xSize(); ++i) 
		for(int j = 0; j < mBounds.ySize()-1; ++j) {
			mBounds(i,j,0) += trans(0);
			mBounds(i,j,1) += trans(1);
			mBounds(i,j,2) += trans(2);
		}
		mCenter += trans;
}

// writeModel
void CMesh::writeModel(const char* aFilename) {
	cout << "Write Model... ";
	std::ofstream aStream(aFilename);
	aStream << "OFF" << std::endl;
	aStream << mNumPoints << " " << mNumPatch << " " << 0 << std::endl;
	// Write vertices
	for (int i = 0; i < mNumPoints; i++) {
		aStream << mPoints[i][0] << " ";   
		aStream << mPoints[i][1] << " ";
		aStream << mPoints[i][2] << std::endl;
	}
	// Write patches
	for (int i = 0; i < mNumPatch; i++) {
		aStream << "3 " << mPatch[i][0] << " ";
		aStream << mPatch[i][1] << " ";
		aStream << mPatch[i][2] << std::endl;
	}
	cout << "ok" << endl;
}

void CMesh::writeAdaptModel(const char* aFilename, const CMesh* adM) {
	cout << "Write Adapt Model... ";
	std::ofstream aStream(aFilename);
	aStream << "OFF" << std::endl;
	aStream << mNumPoints << " " << mNumPatch << " " << 0 << std::endl;
	// Write vertices
	for (int i = 0; i < mNumPoints; i++) {

#if 0
		if(!IsCovered( GetBodyPart(int(mPoints[i][3])) ) ) {
			aStream << mPoints[i][0] << " ";   
			aStream << mPoints[i][1] << " ";
			aStream << mPoints[i][2] << std::endl;
		} else {
			float tmp = mPoints[i][4];
			int j = 5;
			for(;j<mCovered.size()*2+3;) {
				if(IsCovered( GetBodyPart(int(mPoints[i][j])) )) {
					tmp += mPoints[i][++j];
					++j;
				} else j += 2;
			}
			aStream << mPoints[i][0]*(1.0 - tmp) + tmp*adM->mPoints[i][0] << " ";   
			aStream << mPoints[i][1]*(1.0 - tmp) + tmp*adM->mPoints[i][1] << " ";
			aStream << mPoints[i][2]*(1.0 - tmp) + tmp*adM->mPoints[i][2] << std::endl;
		}
#endif    


		float tmp = 0;
		for(int j=3;j<mCovered.size()*2+3;) {
			if(IsCovered( GetBodyPart(int(mPoints[i][j])) )) {
				tmp += mPoints[i][++j];
				++j;
			} else j += 2;
		}

		aStream << mPoints[i][0]*(1.0 - tmp) + tmp*adM->mPoints[i][0] << " ";   
		aStream << mPoints[i][1]*(1.0 - tmp) + tmp*adM->mPoints[i][1] << " ";
		aStream << mPoints[i][2]*(1.0 - tmp) + tmp*adM->mPoints[i][2] << std::endl;
	}
	// Write patches
	for (int i = 0; i < mNumPatch; i++) {
		aStream << "3 " << mPatch[i][0] << " ";
		aStream << mPatch[i][1] << " ";
		aStream << mPatch[i][2] << std::endl;
	}
	cout << "ok" << endl;
}

// writeModel
bool CMesh::adaptOFF(const char* aFilename, float lambda) {
	cout << "Read OFF... ";
	std::ifstream aStream(aFilename);
	if(aStream.is_open()) {
		char buffer[200];
		aStream.getline(buffer,200);
		cout << buffer << endl;
		aStream.getline(buffer,200);
		cout << buffer << endl;
		// Write vertices
		for (int i = 0; i < mNumPoints; i++) {
			for(int j=0; j<3; ++j) {
				float tmp;
				aStream >> tmp;
				mPoints[i][j] *= (1-lambda);
				mPoints[i][j] += lambda*tmp;
			}
		}
		cout << "ok" << endl;
		return true;
	} else return false;
}

// readOFF
bool CMesh::readOFF(const char* aFilename) {
	cout << "Read OFF... ";
	std::ifstream aStream(aFilename);
	if(aStream.is_open()) {

		char buffer[200];
		aStream.getline(buffer,200);
		cout << buffer << endl;
		aStream >> mNumPoints;
		aStream >> mNumPatch;
		aStream >> mNumSmooth;
		mPoints.resize(mNumPoints);
		mPatch.resize(mNumPatch);

		mCenter.setSize(4);
		mCenter[0] = 0;
		mCenter[1] = 0;
		mCenter[2] = 0;
		mCenter[3] = 1.0f;

		mBounds.setSize(1,9,3);
		CVector<float> minV(3,100000);
		CVector<float> maxV(3,-100000);

		// Read vertices
		for (int i = 0; i < mNumPoints; i++) {
			for(int j=0; j<3; ++j) {
				aStream >> mPoints[i][j];
				mCenter[j] += mPoints[i][j];
				if ( mPoints[i][j]<minV(j) ) minV(j) = mPoints[i][j];
				if ( mPoints[i][j]>maxV(j) ) maxV(j) = mPoints[i][j];
			}
		}

		mCenter[0] /= (float)mNumPoints;
		mCenter[1] /= (float)mNumPoints;
		mCenter[2] /= (float)mNumPoints;

		mBounds(0,0,0) = mBounds(0,1,0) = mBounds(0,2,0) = mBounds(0,3,0) = minV(0);
		mBounds(0,4,0) = mBounds(0,5,0) = mBounds(0,6,0) = mBounds(0,7,0) = maxV(0);
		mBounds(0,0,1) = mBounds(0,1,1) = mBounds(0,4,1) = mBounds(0,5,1) = minV(1);
		mBounds(0,2,1) = mBounds(0,3,1) = mBounds(0,6,1) = mBounds(0,7,1) = maxV(1);
		mBounds(0,0,2) = mBounds(0,2,2) = mBounds(0,4,2) = mBounds(0,6,2) = minV(2);
		mBounds(0,1,2) = mBounds(0,3,2) = mBounds(0,5,2) = mBounds(0,7,2) = maxV(2);

		// Read triangles
		for (int i = 0; i < mNumPatch; i++) {
			int dummy;
			aStream >> dummy; 
			for(int j=0; j<3; ++j) {
				aStream >> mPatch[i][j];
			}
		}

		mJointNumber = 0;
		cout << "ok" << endl;
		return true;
	} else return false;
}


void CMesh::writeSkel(const char* aFilename) {
	cout << "Write Skel... ";
	std::ofstream aStream(aFilename);
	aStream << "Skeleton" << std::endl;
	aStream << mJointNumber << std::endl;

	// Write vertices
	for (int i = 1; i <= mJointNumber; i++) {
		aStream << i << " ";   
		aStream << mJoint[i].getDirection() << " ";
		aStream << mJoint[i].getPoint() << " ";
		aStream << mJoint[i].mParent << std::endl;
	}

	cout << "ok" << endl;
}

// rigidMotion
void CMesh::rigidMotion(CVector<CMatrix<float> >& M,CVector<float>& X, bool smooth, bool force) {

	CVector<float> a(4); a(3) = 1.0;
	CVector<float> b(4);
	// Apply motion to points
	if(!smooth || mNumSmooth == 1) {
	  for (int i = 0; i < mNumPoints; i++) {
	    a(0) = mPoints[i](0); a(1) =  mPoints[i](1); a(2) = mPoints[i](2);
	    b = M( int(mPoints[i](3)) )*a;
	    mPoints[i](0)=b(0);mPoints[i](1)=b(1);mPoints[i](2)=b(2);
	  }
	} else {
	  
	  for(int i = 0; i < mNumPoints; i++) {
	    a(0) = mPoints[i](0); a(1) =  mPoints[i](1); a(2) = mPoints[i](2);
	    b = 0;
	    for(int n = 0; n < mNumSmooth; n++)
	      b += M( int(mPoints[i](3 + n * 2)) ) * a * mPoints[i](4 + n * 2);
	    mPoints[i](0)=b(0);mPoints[i](1)=b(1);mPoints[i](2)=b(2);
	    }
	}

	// Apply motion to joints
	for (int i = mJointNumber; i > 0; i--){
	  mJoint(i).rigidMotion(M(mJoint(i).mParent));
	}
	
	if(!smooth || force) {
	  for(int j=0;j<mBounds.xSize();++j) {
	    int jID = (int)mBounds(j,8,0);
	    for(int i=0;i<8;++i) {
	      a(0) = mBounds(j,i,0); a(1) = mBounds(j,i,1); a(2) = mBounds(j,i,2); 
	      b = M(jID)*a;
	      mBounds(j,i,0) = b(0); mBounds(j,i,1) = b(1); mBounds(j,i,2) = b(2);
	    }
	  }
	  
	  mCenter = M(0)*mCenter;
	  
	  // Sum up motion
	  mAccumulatedMotion.mRBM = M(0)*mAccumulatedMotion.mRBM;
	  mCurrentMotion.mRBM = M(0)*mCurrentMotion.mRBM;
	  // Correct pose parameters:
	  CVector<float> T1(6);
	  NRBM::RBM2Twist(T1,mCurrentMotion.mRBM);
	  for (int i=0;i<6;i++)
	    X(i) = T1(i)-mAccumulatedMotion.mPoseParameters(i);
	  mAccumulatedMotion.mPoseParameters += X;
	  mCurrentMotion.mPoseParameters += X;
	}
}

void CMesh::smoothMotionDQ(CVector<CMatrix<float> >& M,CVector<float>& X) {
	std::vector<CVector<float> > vq0(M.size());
	std::vector<CVector<float> > vdq(M.size());

	//transform matrix to quaternion
	for(int i=0; i<M.size(); ++i) {
		vq0[i].setSize(4);
		vdq[i].setSize(4);
		NRBM::RBM2QDQ(M[i],vq0[i],vdq[i]);
	}

	CVector<float> a(4); a(3) = 1.0;
	CVector<float> b(4);
	CVector<float> eq0(4);
	CVector<float> edq(4);
	CMatrix<float> eM(4,4,0);

	// Apply motion to points
	for(int i = 0; i < mNumPoints; i++) {
		a(0) = mPoints[i](0); a(1) =  mPoints[i](1); a(2) = mPoints[i](2);
		eq0 = 0;
		edq = 0;

		CVector<float> q = vq0[ int(mPoints[i](3)) ];
		eq0 = vq0[ int(mPoints[i](3)) ]* mPoints[i](4);
		edq = vdq[ int(mPoints[i](3)) ]* mPoints[i](4);
		for(int n = 1; n < mNumSmooth; n++) {
			if(q*vq0[ int(mPoints[i](3 + n * 2)) ] < 0) {
				eq0 -= vq0[ int(mPoints[i](3 + n * 2)) ] * mPoints[i](4 + n * 2);
				edq -= vdq[ int(mPoints[i](3 + n * 2)) ] * mPoints[i](4 + n * 2);
			} else {
				eq0 += vq0[ int(mPoints[i](3 + n * 2)) ] * mPoints[i](4 + n * 2);
				edq += vdq[ int(mPoints[i](3 + n * 2)) ] * mPoints[i](4 + n * 2);
			}
		}

		float invl = 1.0F/eq0.norm(); 
		eq0 *= invl;
		edq *= invl;

		NRBM::QDQ2RBM(eq0,edq,eM);

		b = eM*a;

		mPoints[i](0)=b(0);mPoints[i](1)=b(1);mPoints[i](2)=b(2);
	}

	// Apply motion to joints
	for (int i = mJointNumber; i > 0; i--) {
		mJoint(i).rigidMotion(M(mJoint(i).mParent));
	}
}


void CMesh::makeSmooth(CMesh* initMesh, bool dual) {
	if(mNumSmooth>1) {

		mPoints = initMesh->mPoints;
		mJoint = initMesh->mJoint;

		for(int i=6; i<mCurrentMotion.mPoseParameters.size(); ++i) {
			while(mCurrentMotion.mPoseParameters[i] < -3.1415926536)
				mCurrentMotion.mPoseParameters[i] += 2.0F*3.1415926536F;
			while(mCurrentMotion.mPoseParameters[i] > 3.1415926536)
				mCurrentMotion.mPoseParameters[i] -= 2.0F*3.1415926536F;
		}

		CVector<float> X(mCurrentMotion.mPoseParameters);
		CVector<CMatrix<float> >M(joints()+1);

		M(0) = mCurrentMotion.mRBM;

		for (int i = 1; i <= joints(); i++) {
			M(i).setSize(4,4); M(i) = 0;
			M(i)(0,0) = 1.0; M(i)(1,1) = 1.0; M(i)(2,2) = 1.0; M(i)(3,3) = 1.0;
		}

		for (int i = joints(); i > 0; i--) {
			CMatrix<float> Mi(4,4);
			joint(i).angleToMatrix(X(5+i),Mi);
			for (int j = 1; j <= joints(); j++) {
				if (influencedBy(j,i)) M(j) = Mi*M(j);
			}
		}

		for (int i = 1; i <= joints(); i++) {
			M(i) = M(0)*M(i);
		}

		if(dual==true)
			smoothMotionDQ(M,X);
		else
			rigidMotion(M,X, true);
	}
}

// angleToMatrix
void CMesh::angleToMatrix(const CMatrix<float>& aRBM, CVector<float>& aJAngles, CVector<CMatrix<float> >& M) {

  // Determine motion of parts behind joints
  
  for (int i = 1; i <= mJointNumber; i++) {
    M(i).setSize(4,4); M(i) = 0;
    M(i)(0,0) = 1.0; M(i)(1,1) = 1.0; M(i)(2,2) = 1.0; M(i)(3,3) = 1.0;
  }
  unsigned int jIds[] = {1,24,2,3,4,5,23,6,7,8,9,10,25,11,12,13,14,15,16,17,18,19,20,21,22};
  
  // Leonid
  for (int ii = mJointNumber; ii > 0; ii--){
    int i = jIds[ii-1];
    CMatrix<float> Mi(4,4);
    mJoint(i).angleToMatrix(aJAngles(i+5),Mi); // i-1
    for (int jj = 1; jj <= mJointNumber; jj++){
      int j = jIds[jj-1];
      if (mInfluencedBy(j,i)) 
	M(j) = Mi*M(j);
    }
  }

  for (int i = 1; i <= mJointNumber; i++)
    M(i) = aRBM*M(i);
  M(0) = aRBM;
}

void CMesh::invAngleToMatrix(const CMatrix<float>& aRBM, CVector<float>& aJAngles, CVector<CMatrix<float> >& M) {
	// Determine motion of parts behind joints
	for (int i = 1; i <= mJointNumber; i++) {
		M(i).setSize(4,4); M(i) = 0;
		M(i)(0,0) = 1.0; M(i)(1,1) = 1.0; M(i)(2,2) = 1.0; M(i)(3,3) = 1.0;
	}
	for (int i = mJointNumber; i > 0; i--) {
		CMatrix<float> Mi(4,4);
		mJoint(i).angleToMatrix(aJAngles(i+5),Mi); // i-1
		for (int j = 1; j <= mJointNumber; j++)
			if (mInfluencedBy(j,i)) M(j) = M(j)*Mi;
	}
	for (int i = 1; i <= mJointNumber; i++)
		M(i) = M(i)*aRBM;
	M(0) = aRBM;
}

void CMesh::twistToMatrix(CVector<float>& aTwist, CVector<CMatrix<float> >& M) {
	NRBM::Twist2RBM(aTwist,M(0));
	angleToMatrix(M(0), aTwist, M);
}

// isParentOf
bool CMesh::isParentOf(int aParentJointID, int aJointID) {
	if (aJointID == 0) return false;
	if (mJoint(aJointID).mParent == aParentJointID) return true;
	return isParentOf(aParentJointID,mJoint(aJointID).mParent);
}

// operator=
void CMesh::operator=(const CMesh& aMesh) {

	mJointNumber = aMesh.mJointNumber;
	mNumPoints=aMesh.mNumPoints;
	mNumPatch=aMesh.mNumPatch;
	mNumSmooth=aMesh.mNumSmooth;

	mPoints = aMesh.mPoints;
	mPatch =  aMesh.mPatch;

	mNoOfBodyParts = aMesh.mNoOfBodyParts;
	mBoundJoints = aMesh.mBoundJoints;


	mBounds=aMesh.mBounds;
	mCenter=aMesh.mCenter;
	mJointMap=aMesh.mJointMap;
	mNeighbor=aMesh.mNeighbor;
	mEndJoint=aMesh.mEndJoint;

	mCovered=aMesh.mCovered;
	mExtremity=aMesh.mExtremity;

	mJoint = aMesh.mJoint;
	mInfluencedBy = aMesh.mInfluencedBy;

	mAccumulatedMotion = aMesh.mAccumulatedMotion;
	mCurrentMotion = aMesh.mCurrentMotion;
}

void CMesh::projectToImage(CMatrix<float>& aImage, CMatrix<float>& P, int aLineValue)
{

	int I[4];
	CVector<float> aPoint(3);
	int ax,ay,bx,by,cx,cy;
	int Size = (*this).GetMeshSize();
	int j;

	int det;

	float X,Y,Z;

	//cout << Size << endl;

	for(int i=0;i<Size;i++)
	{
		(*this).GetPatch(i, I[0],I[1],I[2]);
		j=0;
		// Test if patch is visible ...
		X= mPoints[ I[0] ](0);Y=mPoints[ I[0] ](1);Z=mPoints[ I[0] ](2);
		(*this).projectPoint(P, X,Y,Z,ax,ay   );
		X= mPoints[ I[1] ](0);Y=mPoints[ I[1] ](1);Z=mPoints[ I[1] ](2);
		(*this).projectPoint(P, X,Y,Z,bx,by   );
		X= mPoints[ I[2] ](0);Y=mPoints[ I[2] ](1);Z=mPoints[ I[2] ](2);
		(*this).projectPoint(P, X,Y,Z,cx,cy   );

		det=ax*(by-cy) - bx*(ay-cy) + cx*(ay-by);

		if(det<0.001)
		{
			for(j=0;j<3;j++)
			{
				(*this).projectPoint(P, mPoints[ I[j] ](0),mPoints[ I[j] ](1),mPoints[ I[j] ](2),bx,by);


				if (ax >= 0 && ay >= 0 && ax < aImage.xSize() && ay < aImage.ySize())
					aImage.drawLine(ax,ay,bx,by,(float)aLineValue);
				ax=bx;
				ay=by;
			}
			j=0;
			(*this).projectPoint(P, mPoints[ I[j] ](0),mPoints[ I[j] ](1),mPoints[ I[j] ](2),bx,by);
			if (ax >= 0 && ay >= 0 && ax < aImage.xSize() && ay < aImage.ySize())
				aImage.drawLine(ax,ay,bx,by,(float)aLineValue);
		}
	}
}

void CMesh::projectToImage(CMatrix<float>& aImage, CMatrix<float>& P, int aLineVal, int aNodeVal) {
	projectToImage(aImage, P,  aLineVal);
	projectPointsToImage(aImage,P,aNodeVal);
}

void CMesh::projectToImage(CTensor<float>& aImage, CMatrix<float>& P, int aLineR, int aLineG, int aLineB) {
	CMatrix<float> aMesh(aImage.xSize(),aImage.ySize(),0);
	projectToImage(aMesh,P,255);
	int aSize = aMesh.size();
	int a2Size = 2*aSize;
	for (int i = 0; i < aSize; i++)
		if (aMesh.data()[i] == 255) {
			aImage.data()[i] = (float)aLineR;
			aImage.data()[i+aSize] = (float)aLineG;
			aImage.data()[i+a2Size] = (float)aLineB;
		}
#if 0
		int mx, my;
		projectPoint(P, mCenter[0], mCenter[1], mCenter[2],mx,my);
		int crossSize = 10;
		int x1 = NMath::max(mx-crossSize,0); int x2 = NMath::min(mx+crossSize,aImage.xSize()-1);
		int y1 = NMath::max(my-crossSize,0); int y2 = NMath::min(my+crossSize,aImage.ySize()-1);
		for(int x = x1; x<=x2; ++x) {
			aImage(x,my,0) = 40;
			aImage(x,my,1) = 240;
			aImage(x,my,2) = 40;
		}
		for(int y = y1; y<=y2; ++y) {
			aImage(mx,y,0) = 40;
			aImage(mx,y,1) = 240;
			aImage(mx,y,2) = 40;
		}
#endif
}

void CMesh::projectToImage(CTensor<float>& aImage, CMatrix<float>& P, int aLineR, int aLineG, int aLineB, int aNodeR, int aNodeG, int aNodeB) {
	CMatrix<float> aMesh(aImage.xSize(),aImage.ySize(),0);
	projectToImage(aMesh,P,255,128);
	int aSize = aMesh.size();
	int a2Size = 2*aSize;
	for (int i = 0; i < aSize; i++)
		if (aMesh.data()[i] == 255) {
			aImage.data()[i] = (float)aLineR;
			aImage.data()[i+aSize] = (float)aLineG;
			aImage.data()[i+a2Size] = (float)aLineB;
		}
		else if (aMesh.data()[i] == 128) {
			aImage.data()[i] = (float)aNodeR;
			aImage.data()[i+aSize] = (float)aNodeG;
			aImage.data()[i+a2Size] = (float)aNodeB;
		}
}

// projectPointsToImage
void CMesh::projectPointsToImage(CMatrix<float>& aImage, CMatrix<float>& P, int aNodeVal) {
	int ax,ay;
	CMatrix<float> aMesh(aImage.xSize(),aImage.ySize(),0);
	projectToImage(aMesh, P,  aNodeVal);

	int Size = (*this).GetPointSize();
	for(int i=0; i< Size;i++)
	{
		projectPoint(P,mPoints[i](0),mPoints[i](1),mPoints[i](2),ax,ay);
		if (ax >= 0 && ay >= 0 && ax < aImage.xSize() && ay < aImage.ySize())
			if(aMesh(ax,ay)==aNodeVal)
				aImage(ax,ay) = (float)aNodeVal;
	}
}


void CMesh::projectPointsToImage(CTensor<float>& a3DCoords, CMatrix<float>& P) {
	int ax,ay;
	int aNodeVal=39;

	CMatrix<float> aMesh(a3DCoords.xSize(),a3DCoords.ySize(),0);
	projectToImage(aMesh, P,  aNodeVal);

	int Size = GetPointSize();
	for (int i = 0; i < Size; i++) 
	{
		projectPoint(P,mPoints[i](0),mPoints[i](1),mPoints[i](2),ax,ay);
		if (ax >= 0 && ay >= 0 && ax < a3DCoords.xSize() && ay < a3DCoords.ySize())
			if(aMesh(ax,ay)==aNodeVal)
			{
				a3DCoords(ax,ay,0) = mPoints[i](0);
				a3DCoords(ax,ay,1) = mPoints[i](1);
				a3DCoords(ax,ay,2) = mPoints[i](2);
				a3DCoords(ax,ay,3) = mPoints[i](3);
			}
	}

}

void CMesh::projectToImageJ(CMatrix<float>& aImage, CMatrix<float>& P, int aLineValue, int aJoint, int xoffset, int yoffset)
{

	int I[4];
	CVector<float> aPoint(3);
	int ax,ay,bx,by,cx,cy;
	int Size = (*this).GetMeshSize();
	int j;

	int det;

	float X,Y,Z;

	for(int i=0;i<Size;i++)
	{
		(*this).GetPatch(i, I[0],I[1],I[2]);
		j=0;
		// Test if patch is visible ...
		X= mPoints[ I[0] ](0);Y=mPoints[ I[0] ](1);Z=mPoints[ I[0] ](2);
		(*this).projectPoint(P, X,Y,Z,ax,ay   );
		ax -= xoffset; ay -= yoffset;
		X= mPoints[ I[1] ](0);Y=mPoints[ I[1] ](1);Z=mPoints[ I[1] ](2);
		(*this).projectPoint(P, X,Y,Z,bx,by   );
		bx -= xoffset; by -= yoffset;
		X= mPoints[ I[2] ](0);Y=mPoints[ I[2] ](1);Z=mPoints[ I[2] ](2);
		(*this).projectPoint(P, X,Y,Z,cx,cy   );
		cx -= xoffset; cy -= yoffset;

		det=ax*(by-cy) - bx*(ay-cy) + cx*(ay-by);

		{
			if(  ((*this).GetJointID(I[0])==aJoint) || ((*this).GetJointID(I[1])==aJoint) || ((*this).GetJointID(I[2])==aJoint))
			{
				for(j=0;j<3;j++)
				{
					(*this).projectPoint(P, mPoints[ I[j] ](0),mPoints[ I[j] ](1),mPoints[ I[j] ](2),bx,by);
					bx -= xoffset; by -= yoffset;

					if (ax >= 0 && ay >= 0 && ax < aImage.xSize() && ay < aImage.ySize())
						aImage.drawLine(ax,ay,bx,by,(float)aLineValue);
					ax=bx;
					ay=by;
				}
				j=0;
				(*this).projectPoint(P, mPoints[ I[j] ](0),mPoints[ I[j] ](1),mPoints[ I[j] ](2),bx,by);
				bx -= xoffset; by -= yoffset;
				if (ax >= 0 && ay >= 0 && ax < aImage.xSize() && ay < aImage.ySize())
					aImage.drawLine(ax,ay,bx,by,(float)aLineValue);
			}
		}
	}
}

// projectSurface
void CMesh::projectSurface(CMatrix<float>& aImage, CMatrix<float>& P, int xoffset, int yoffset) {
	//aImage = 0;

	CMatrix<float> ImageP=aImage;

	for (int k = 0; k <= (*this).joints(); k++)
	{ 
		ImageP = 0;
		projectToImageJ(ImageP,P,1,k, xoffset, yoffset);

		for(int i=0;i<ImageP.xSize();i++)
			ImageP(i,0)=0;
		for(int i=0;i<ImageP.ySize();i++)
			ImageP(0,i)=0;
		int aSize = ImageP.size();
		for (int i = 0; i < aSize; i++)
			if (ImageP.data()[i] == 0) {
				int y = i/ImageP.xSize();
				int x = i-ImageP.xSize()*y;
				ImageP.connectedComponent(x,y);
				break;
			}
			// Fill in ...
			for (int i = 0; i < aSize; i++)
			{
				if (ImageP.data()[i] == 0) aImage.data()[i] = 255;
			}
	}
}

void CMesh::readShapeSpaceEigens(const double* eigenVectorsIn, int numEigenVectors, int nPoints){
  
  eigenVectors.resize(numEigenVectors);
  int dimD = 3;
  
  for (int i0 = 0; i0 < numEigenVectors; i0++){
    eigenVectors[i0].setSize(nPoints,dimD);
  }
  int n = 0;
  for(int col = 0; col < dimD; col++)
    for(int row = 0; row < nPoints; row++)
      for (int i0 = 0; i0 < numEigenVectors; i0++){
	eigenVectors[i0](row, col) = eigenVectorsIn[n];
	n++;
      }
}

void CMesh::readShapeSpaceEigens(std::string fileName, int numEigenVectors)
{
	string path = fileName;		
	ifstream loadFile;
	loadFile.open(path.c_str());
	//reading only first four columns here
	//char fullLine[10000], c[numEigenVectors][500];
	char* fullLine = new char[10000];
	char** c = new char* [numEigenVectors];
	for(int i = 0; i<numEigenVectors; i++)
		c[i] = new char[500];
	eigenVectors.resize(numEigenVectors);
	for (int i0 = 0; i0 < numEigenVectors; i0++)
	{
		eigenVectors[i0].setSize(13249,3);
		eigenVectors[i0] = 0;
	}
	int row = 0, col = 0;
	unsigned totalLines = 1;
	int noPts = 0;
	while(!loadFile.getline(fullLine, 10000).eof())
	{
		if(fullLine[0] == '#') 
		{	//first start line with comments
			int totalRws = 0, dimD = 0;
			sscanf(fullLine,"# %d ( %d x %d )", &totalRws, &noPts, &dimD);
			continue;
		}
		sscanf(fullLine,"%s%s%s%s%s%*s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",
			&c[0], &c[1], &c[2], &c[3], 
			&c[4], &c[5], &c[6], &c[7],&c[8], &c[9],
			&c[10], &c[11], &c[12], &c[13], 
			&c[14], &c[15], &c[16], &c[17],&c[18], &c[19]);	
		for (int i0 = 0; i0 < numEigenVectors; i0++)
		{	

			eigenVectors[i0](row, col) = atof(c[i0]);
		}
		row++; 
		if(row >= noPts)
		{
			row = 0; col++;
			if(col == 3) col = 0;
		}
		totalLines++;
	}
	loadFile.close();		
}

int CMesh::shapeChangesToMesh( CVector<float> shapeParams ){
    
  size_t noPts = mPoints.size();	
  for(unsigned int i1 = 0; i1 < mPoints.size(); i1++){
    for(unsigned int i2 = 0; i2 < 3; i2++){
      for(int i0 = 0; i0 < shapeParams.size(); i0++){					
	double tmp = shapeParams[i0];
	(mPoints[i1])[i2] += (float)shapeParams[i0]*(float)SKEL_SCALE_FACTOR*(float)eigenVectors[i0](i1,i2);
      }
    }
  }	
  //cout<<"Done!";
  return true;
}

int CMesh::updateJntPos()
{
  // cout << "\nUpdating Joints.. ";	
  CMatrix<float> newJntPos;newJntPos.setSize(mJointNumber, 3); //3D joints
	CVector<float> newJnt; newJnt.setSize(3);
	CMatrix<float> tmpMatrix; tmpMatrix.setSize(weightMatrix.xSize(), weightMatrix.ySize()); tmpMatrix = 0;
	CMatrix<float> joints0; joints0.setSize(mJointNumber, mJointNumber*3);
	CVector<float> minEle; minEle.setSize(mNumPoints);
	CVector<float> singleWts; singleWts.setSize(mNumPoints); singleWts = 0;

	for(int i0 = 0; i0 < mJointNumber; i0++)
	  {
	    for(int i1 = 0; i1 < mJointNumber; i1++)
	      {
		if(
		   
		   (i0 == 2 && i1 ==0) || (i0 == 3 && i1 == 2) || (i0 == 4 && i1 == 3)||
		   (i0 == 6 && i1 == 0) || (i0 == 7 && i1 == 6) || (i0 == 8 && i1 == 7) ||
		   (i0 == 10 && i1 == 0) || (i0 == 11 && i1 == 10) || (i0 == 14 && i1 == 10) ||
		   (i0 == 15 && i1 == 14) || (i0 == 16 && i1 == 15) || (i0 == 19 && i1 == 10) ||
		   (i0 == 20 && i1 == 19) || (i0 == 21 && i1 == 20) 
		   
		   )
		  {
		    
		    bool anyMinIsNotZero = minMLab(weightMatrix, i0, i1, minEle);
		    
		    if(anyMinIsNotZero)
		      {	//printf("\nNon-zero min wt bw joints %d and %d", i0+1, i1+1);
			componet_wise_mul_with_pnts(minEle, tmpMatrix);
			double sumMinEle = sumTheVector(minEle); assert(sumMinEle > 0);			
			CVector<float> t1; t1.setSize(tmpMatrix.ySize());
			
			sumTheMatrix(tmpMatrix, t1);				
			
			newJnt(0) = t1(0)/(float)sumMinEle;
			newJnt(1) = t1(1)/(float)sumMinEle;
			newJnt(2) = t1(2)/(float)sumMinEle;
			
			joints0(i0, 3*i1 + 0) = newJnt(0);
			joints0(i0, 3*i1 + 1) = newJnt(1);
			joints0(i0, 3*i1 + 2) = newJnt(2);
			
		      }
		  }
	      }
	  }	
	copyMatColToVector(weightMatrix, 0, singleWts);
	tmpMatrix = 0;
	componet_wise_mul_with_pnts(singleWts, tmpMatrix);
	CVector<float> t1; t1.setSize(tmpMatrix.ySize());
	sumTheMatrix(tmpMatrix, t1);
	double sumRootElse = sumTheVector(singleWts); assert(sumRootElse > 0);
	if(sumRootElse <= 0) sumRootElse = 1;
	newJntPos(0, 0) = t1(0)/(float)sumRootElse;
	newJntPos(0, 1) = t1(1)/(float)sumRootElse;
	newJntPos(0, 2) = t1(2)/(float)sumRootElse;

	int parentMap[50]; // jNo=(jId -1) (to)---> parent jId
	for(int i0 = 0; i0 < mJointNumber; i0++)
	  {
	    int jId = i0 + 1;
	    parentMap[i0] = mJoint(jId).mParent;
	  }
	
	parentMap[2] = 1;
	parentMap[3] = 3;
	parentMap[5] = 5; 
	parentMap[6] = 1; 
	parentMap[7] = 7;
	parentMap[8] = 8;
	parentMap[9] = 9; 
	parentMap[10] = 1; 
	parentMap[14] = 11;
	parentMap[19] = 11;
	parentMap[22] = 22; 
	parentMap[23] = 23; 
	parentMap[24] = 24; 
	parentMap[25] = 25;
	
	std::set<int> nonRoots; //has jId's
	
	nonRoots.insert(3);nonRoots.insert(4);
	nonRoots.insert(5);nonRoots.insert(7);
	nonRoots.insert(8);nonRoots.insert(9);
	nonRoots.insert(11);nonRoots.insert(12);
	nonRoots.insert(15);nonRoots.insert(16);
	nonRoots.insert(17);nonRoots.insert(20);
	nonRoots.insert(21);nonRoots.insert(22);
	nonRoots.insert(23);nonRoots.insert(24);
	nonRoots.insert(25);

	float tmp[3];
	for(int i0 = 0; i0 < mJointNumber; i0++)
	  {
	    int jId = i0 + 1;		
	    if(nonRoots.count(jId) != 0)
	      {	//non root joints
		
		tmp[0] = newJntPos(i0, 0) = joints0( i0, (parentMap[i0] - 1)*3 + 0);
		tmp[1] = newJntPos(i0, 1) = joints0( i0, (parentMap[i0] - 1)*3 + 1);
		tmp[2] = newJntPos(i0, 2) = joints0( i0, (parentMap[i0] - 1)*3 + 2);	
	      }
	  }
	copyJointPos(1, 2, newJntPos);
	copyJointPos(5, 6, newJntPos);
	copyJointPos(9, 10, newJntPos);	
	copyJointPos(12, 14, newJntPos);
	copyJointPos(13, 14, newJntPos);
	copyJointPos(17, 19, newJntPos);
	copyJointPos(18, 19, newJntPos);
	copyJointPos(22, 6, newJntPos);
	copyJointPos(23, 2, newJntPos);
	copyJointPos(24, 10, newJntPos);

	for(int i0 = 0; i0 < mJointNumber; i0++)
	  {
	    unsigned int jId = i0 + 1;
	    CVector<float> aDir = mJoint(jId).getDirection();		
	    CVector<float> jPos; jPos.setSize(3);
	    
	    jPos(0) = newJntPos(i0,0);
	    jPos(1) = newJntPos(i0,1);
	    jPos(2) = newJntPos(i0,2);
	    mJoint(jId).set(aDir, jPos);
	  }
	return 1;
}

bool CMesh::minMLab( CMatrix<float> weightMatrix, int i0, int i1, CVector<float> &minEle )
{
	bool isNotZero = false;
	for(int i2 = 0; i2 < weightMatrix.xSize(); i2++)
	{		
		float w1 = weightMatrix(i2,i0);
		float w2 = weightMatrix(i2, i1);
		if(weightMatrix(i2,i0) < weightMatrix(i2, i1))
			minEle(i2) = weightMatrix(i2,i0);
		else
			minEle(i2) = weightMatrix(i2, i1);
		if(minEle(i2) > 0) isNotZero = true;
	}	
	return isNotZero;
}

void CMesh::componet_wise_mul_with_pnts( CVector<float> minEle, CMatrix<float> &tmpMatrix )
{
	assert(mPoints.size() == minEle.size());
	bool allZero = true;
	for(int i0 = 0; i0 < mPoints.size(); i0++)
	{
		if(minEle(i0) > 0) 
		{
			double m_ele = minEle(i0);
			allZero = false;
		}
		float tmp[3];
		for(int i1 = 0; i1 < 3; i1++)
		{			
			tmp[i1] = (mPoints[i0])(i1);
			tmpMatrix(i0, i1) = (mPoints[i0])(i1) * minEle(i0);
		}
	}
	if(allZero == true)
		printf("\nSomeing the matter, all weights were zero.");
}

double CMesh::sumTheVector( CVector<float> minEle )
{
	double sum = 0.;
	for(int i0 = 0; i0 < minEle.size(); i0++)
		sum += minEle(i0);
	return sum;
}

void CMesh::sumTheMatrix( CMatrix<float> tmpMatrix, CVector<float> & t1 )
{
	t1 = 0.;
	for(int i0 = 0; i0 < tmpMatrix.xSize(); i0++)
	{
		for(int i1 = 0; i1 < tmpMatrix.ySize(); i1++)
		{
			t1(i1) += tmpMatrix(i0, i1);
		}
	}
}

bool CMesh::copyMatColToVector( CMatrix<float> weightMatrix, int col, CVector<float> &singleWts )
{
	bool isNotZero = false;
	double total = 0;
	for(int i0 = 0; i0 < weightMatrix.xSize(); i0++)
	{
		singleWts(i0) = weightMatrix(i0, col);
		if(singleWts(i0) > 0) isNotZero = true;
		total += weightMatrix(i0, col);
	}
	return isNotZero;
}

void CMesh::copyJointPos( int to, int from, CMatrix<float> &newJntPos )
{
	newJntPos(to, 0) = newJntPos(from, 0); 
	newJntPos(to, 1) = newJntPos(from, 1); 
	newJntPos(to, 2) = newJntPos(from, 2);
}

void CMesh::printPoints( std::string fname )
{
	NShow::mLogFile.open( (NShow::mResultDir + fname).c_str(), std::ios_base::app );
	for(unsigned int i0 = 0; i0 < mPoints.size(); i0++)
	{
		NShow::mLogFile << mPoints[i0](0) <<" " << 
			mPoints[i0](1)<< " " <<
			mPoints[i0](2)<<"\n";
	}

	NShow::mLogFile.close();
}

void CMesh::findNewAxesOrient( CVector<double> &axisLeftArm, CVector<double> &axisRightArm, CMatrix<float> newJntPos )
{
	CVector<double> upperArmL, upperArmR, lowerArmL, lowerArmR;
	upperArmL.setSize(3); lowerArmL.setSize(3);
	upperArmR.setSize(3); lowerArmR.setSize(3);
#define WRISTPTS_R 307
#define WRISTPTS_L 424
	unsigned int rightWristPts[WRISTPTS_R] = { 386, 382, 380, 379, 375, 373, 372, 371, 370, 368, 366, 346, 364, 363, 369, 360, 357, 362, 353, 352, 356, 347, 343, 340, 339, 338, 336, 342, 335, 334, 333, 330, 329, 324, 321, 318, 317, 314, 313, 312, 319, 310, 305, 304, 303, 302, 300, 299, 296, 298, 295, 294, 290, 284, 289, 283, 281, 280, 278, 277, 276, 275, 274, 269, 268, 266, 265, 264, 261, 259, 258, 257, 255, 254, 251, 250, 249, 248, 252, 245, 244, 241, 240, 239, 235, 237, 233, 232, 231, 228, 230, 224, 223, 222, 220, 219, 226, 217, 216, 213, 212, 214, 211, 210, 209, 208, 207, 206, 205, 203, 199, 198, 197, 195, 194, 193, 192, 191, 202, 190, 189, 187, 186, 185, 183, 182, 188, 181, 180, 179, 176, 184, 175, 173, 172, 171, 170, 169, 168, 174, 165, 164, 163, 162, 166, 160, 159, 158, 167, 157, 155, 154, 152, 153, 151, 150, 149, 161, 148, 147, 146, 156, 145, 144, 143, 142, 141, 139, 138, 140, 136, 137, 135, 133, 132, 130, 134, 129, 128, 127, 126, 125, 124, 123, 122, 121, 119, 118, 120, 117, 131, 116, 115, 114, 113, 112, 111, 109, 108, 107, 106, 105, 103, 102, 101, 104, 100, 99, 110, 98, 97, 95, 96, 94, 93, 92, 91, 90, 89, 88, 87, 86, 85, 84, 83, 81, 80, 79, 78, 82, 77, 76, 74, 73, 71, 72, 70, 75, 69, 67, 66, 68, 65, 64, 60, 59, 62, 58, 57, 56, 63, 55, 61, 54, 52, 51, 50, 53, 49, 48, 47, 46, 43, 42, 40, 45, 38, 41, 37, 35, 34, 44, 31, 36, 39, 30, 33, 29, 32, 28, 27, 26, 25, 24, 22, 21, 20, 19, 17, 23, 12, 13, 18, 11, 16, 10, 15, 9, 5, 7, 8, 6, 4, 14, 1, 3, 2 };

	unsigned int leftWristPts[WRISTPTS_L] = { 11512, 11514, 11513, 11506, 11508, 11505, 11510, 11501, 11511, 11500, 11507, 11497, 11509, 11492, 11491, 11503, 11504, 11496, 11499, 11484, 11490, 11487, 11498, 11489, 11493, 11477, 11502, 11495, 11494, 11480, 11469, 11483, 11488, 11476, 11465, 11464, 11463, 11485, 11468, 11459, 11475, 11474, 11478, 11458, 11455, 11486, 11481, 11467, 11462, 11448, 11456, 11479, 11446, 11461, 11450, 11457, 11447, 11442, 11482, 11471, 11473, 11440, 11454, 11466, 11435, 11472, 11452, 11470, 11428, 11449, 11451, 11424, 11443, 11439, 11427, 11432, 11422, 11418, 11436, 11425, 11426, 11415, 11438, 11444, 11429, 11412, 11460, 11409, 11445, 11431, 11405, 11421, 11403, 11402, 11453, 11430, 11400, 11411, 11434, 11441, 11396, 11420, 11417, 11394, 11416, 11392, 11419, 11391, 11408, 11389, 11388, 11387, 11413, 11437, 11382, 11398, 11381, 11393, 11379, 11395, 11378, 11376, 11375, 11410, 11374, 11414, 11433, 11423, 11373, 11383, 11386, 11368, 11377, 11380, 11401, 11364, 11397, 11407, 11390, 11385, 11365, 11366, 11370, 11346, 11404, 11344, 11406, 11343, 11339, 11369, 11336, 11350, 11353, 11347, 11341, 11332, 11354, 11399, 11345, 11334, 11349, 11326, 11355, 11371, 11328, 11324, 11384, 11361, 11372, 11333, 11319, 11338, 11363, 11356, 11362, 11314, 11313, 11340, 11357, 11306, 11331, 11309, 11305, 11367, 11360, 11335, 11300, 11351, 11359, 11358, 11301, 11342, 11352, 11299, 11311, 11295, 11348, 11323, 11294, 11327, 11308, 11316, 11281, 11320, 11310, 11318, 11297, 11271, 11315, 11321, 11272, 11280, 11263, 11262, 11296, 11291, 11261, 11330, 11337, 11329, 11254, 11322, 11268, 11307, 11278, 11302, 11286, 11287, 11240, 11284, 11282, 11220, 11234, 11290, 11244, 11224, 11269, 11298, 11223, 11249, 11277, 11270, 11218, 11245, 11216, 11279, 11233, 11248, 11222, 11207, 11241, 11236, 11206, 11205, 11260, 11204, 11199, 11229, 11253, 11197, 11195, 11274, 11258, 11273, 11186, 11237, 11208, 11252, 11250, 11219, 11178, 11264, 11257, 11173, 11185, 11184, 11169, 11239, 11198, 11183, 11196, 11230, 11179, 11242, 11201, 11161, 11228, 11167, 11213, 11176, 11174, 11166, 11182, 11158, 11165, 11192, 11157, 11162, 11146, 11212, 11210, 11144, 11164, 11150, 11163, 11140, 11147, 11203, 11180, 11135, 11194, 11181, 11130, 11143, 11154, 11125, 11137, 11129, 11145, 11133, 11121, 11118, 11123, 11134, 11115, 11188, 11168, 11139, 11107, 11127, 11106, 11155, 11114, 11109, 11111, 11096, 11159, 11105, 11097, 11128, 11085, 11122, 11083, 11112, 11084, 11124, 11061, 10997, 11142, 11069, 11100, 11067, 11081, 11098, 11116, 11070, 11056, 11091, 11054, 11104, 11095, 11051, 11080, 11076, 11066, 11071, 11065, 11045, 11088, 11049, 11043, 11055, 11089, 11040, 11039, 11092, 11029, 11041, 11032, 11028, 11048, 11046, 11050, 11053, 11026, 11060, 11010, 11020, 11017, 11031, 11011, 11006, 11022, 11025, 10998, 11018, 11016, 11000, 10996, 11003, 10991, 10987, 11019, 10995, 10986, 10985, 11021, 11015, 10977, 10974, 10982, 10999, 10984, 10943, 10939, 10965, 10961, 10970, 10950, 10951, 10933, 10976, 10949, 10928 };
	CVector<double> leftWrist, rightWrist;
	leftWrist.setSize(3); rightWrist.setSize(3);
	leftWrist = rightWrist = 0;

	for(unsigned int i0 = 0; i0 < WRISTPTS_R; i0++)
	{
		unsigned int rIdx = rightWristPts[i0];
		for(unsigned int i1 = 0; i1 < 3; i1++)
		{		
			rightWrist(i1) += (mPoints[rIdx])(i1);
		}		
	}

	for (unsigned int i0 = 0; i0 < WRISTPTS_L; i0++)
	{
		unsigned int lIdx = leftWristPts[i0];
		for (unsigned int i1 = 0; i1 < 3; i1++)
		{
			leftWrist(i1) += (mPoints[lIdx])(i1);
		}
	}

	for(unsigned int i1 = 0; i1 < 3; i1++)
	{
		leftWrist(i1) = leftWrist(i1)/WRISTPTS_L;
		rightWrist(i1) = rightWrist(i1)/ WRISTPTS_R;
	}		

	for(unsigned i0 = 0; i0 < 3; i0++)
	{
		upperArmR(i0) = newJntPos(7,i0) - newJntPos(6,i0);
		lowerArmR(i0) = rightWrist(i0) - newJntPos(7,i0);
		upperArmL(i0) = newJntPos(11,i0) - newJntPos(10,i0);
		lowerArmL(i0) = leftWrist(i0) - newJntPos(11,i0);
	}
	double normUpArmR = upperArmR.norm(); 
	double normLowArmR = lowerArmR.norm();
	double normUpArmL = upperArmL.norm();
	double normLowArmL = lowerArmL.norm();
	for(unsigned i0 = 0; i0 < 3; i0++)
	{
		upperArmR(i0) /= normUpArmR;
		lowerArmR(i0) /= normLowArmR;
		upperArmL(i0) /= normUpArmL;
		lowerArmL(i0) /= normLowArmL;
	}
	//////////////////////////////////////////////////////////////////////////
	axisLeftArm = upperArmL / lowerArmL;
	axisRightArm = upperArmR / lowerArmR ;

	double normAxisLeftArm = axisLeftArm.norm();
	double normAxisRightArm = axisRightArm.norm();

	if (normAxisLeftArm == 0 || axisRightArm == 0)
	{
		printf("\n!!!!--- WORRY WORRY AXIS 0 ---!!!");
		assert(0);
	}
	for(unsigned i0 = 0; i0 < 3; i0++)
	{
		axisLeftArm(i0) /= normAxisLeftArm;
		axisRightArm(i0) /= normAxisRightArm;
	}
}

void CMesh::writeMeshDat( std::string fname )
{
	NShow::mLogFile.open( fname.c_str() );
	char buffer[10000];	
	//sprintf(buffer,"%d %d %d %d\n", mPoints.size(), mPatch.size(), joints(), mNoOfBodyParts );
	NShow::mLogFile << buffer;
	for (int i0 = 0; i0 < mNumPoints; i0++ )
	{
		sprintf(buffer,"%f %f %f", mPoints[i0](0), mPoints[i0](1), mPoints[i0](2));
		NShow::mLogFile << buffer;
		getWeightsinBuffer(buffer, i0);
		NShow::mLogFile << buffer <<"\n";
	}
	for (unsigned int i0 = 0; i0 < mPatch.size(); i0++ )
	{
		sprintf(buffer,"%d %d %d\n", mPatch[i0](0), mPatch[i0](1), mPatch[i0](2));
		NShow::mLogFile << buffer;
	}
	for(int k0 = 0; k0 < mJointNumber; k0++)
	{
		int jId = k0 + 1;
		CVector<float> &t_jPos = mJoint(jId).getPoint();
		CVector<float> &t_jDir = mJoint(jId).getDirection();
		sprintf(buffer, "%d %f %f %f %f %f %f %d\n", jId,
			t_jDir(0), t_jDir(1), t_jDir(2), t_jPos(0), t_jPos(1), t_jPos(2),
			mJoint(jId).mParent);
		NShow::mLogFile << buffer;
	}
	NShow::mLogFile.close();
	cout<<"\bDone Writing File!"<<fname;
}

void CMesh::getWeightsinBuffer( char *buffer, int ptNo )
{
	std::vector< std::pair<double, int> > wts;
	std::pair<double, int> tmp;
	string wtStr = "";
	for(int i1 = joints() - 1; i1 >= 0; i1--)
	{
		if(mBoundJoints(i1+1) == true)		
			wts.push_back(make_pair(weightMatrix(ptNo, i1),i1));		
	}

	for ( int i = 0; i < wts.size(); ++i )
	{
		for ( int j = 1; j < wts.size() - i; ++j )
		{
			if ( wts[j-1].first < wts[j].first ) 
			{
				tmp = wts[j-1];
				wts[j-1] = wts[j];
				wts[j] = tmp;
			}			
		}
	}
	
	for(int i0 = 0; i0 < wts.size(); i0++)
	{
		sprintf(buffer, " %d %f", wts[i0].second + 1, wts[i0].first);
		wtStr+= buffer;
	}
	strcpy(buffer, wtStr.c_str());	
}
