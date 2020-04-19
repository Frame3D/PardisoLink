#include <algorithm>
#include "LTemplate.h"

extern "C" void pardiso (void *, mint *, mint *, mint *, mint *, mint *, mreal *, mint *, mint *, mint *, mint *, mint *, mint *, mreal *, mreal *, mint *);
class Pardiso
{
	mint i, n, nnz, nrhs, mtype, maxfct, mnum, phase, error, msglvl;
    mint *ia, *ja, *perm;
	mreal* a;
	mint iparm[64];
	void* pt[64];
	mreal ddum;
	mint idum;
	mreal timestamp;
	mint initialized = 0;
	mint symfactorized = 0;
	mint numfactorized = 0;
    mint checkmatrix = 1;
	
	public:
	;
	Pardiso()
	{
        for ( i = 0; i < 64; i++ )
        {
            iparm[i] = 0;
        }
        for ( i = 0; i < 64; i++ )
        {
            pt[i] = 0;
        }
	};
	~Pardiso()
	{
		if(initialized)
		{
			phase = -1;
			pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
		}
        if(a)
        {
            free(a);
        }
        if(ia)
        {
            free(ia);
        }
        if(ja)
        {
            free(ja);
        }
        if(perm)
        {
            free(perm);
        }
	};
	mma::IntTensorRef RowPointers()
	{
        if(ia)
        {
            return mma::makeVector<mint>(n+1,ia);
        }
	}
	mma::IntTensorRef ColumnIndices()
	{
        if(ja)
        {
            return mma::makeVector<mint>(nnz,ja);
        }
	}
	mma::RealTensorRef NonzeroValues()
	{
        if(a)
        {
            return mma::makeVector<mreal>(nnz,a);
        }
	}
	mreal TimeStamp()
	{
		return timestamp;
	}
	mma::IntTensorRef IntegerParameters()
	{
		return mma::makeVector<mint>(64,iparm);
	}
	mint NumberOfNonzeros()
	{
		return nnz;
	}
	mint Length()
	{
		return n;
	}
	mint Error()
	{
		return error;
	}
	mint MatrixType()
	{
		return mtype;
	}
    mint CheckMatrixQ()
    {
        return checkmatrix;
    }
	mma::IntTensorRef Permutation()
	{
        if(perm)
        {
            return mma::makeVector<mint>(n,perm);
        }
	}
	void Init(mma::IntTensorRef ia0, mma::IntTensorRef ja0, mma::RealTensorRef a0, mint mtype0)
	{
        if(!initialized)
        {
            mtype = mtype0;

            //pardisoinit(pt,&mtype,iparm);

            nrhs = 1;
            timestamp = 0.;
            n = ia0.length()-1;
            nnz = ia0[ia0.length()-1]-1;
            ia = (mint*) malloc(ia0.length()*sizeof(mint));
            if(!ia)
            {
                mma::print("Allocation of ia failed.");
                error -100001;
            }
            std::copy(ia0.data(),ia0.data()+ia0.length(),ia);
            ja = (mint*) malloc(ja0.length()*sizeof(mint));
            if(!ja)
            {
                mma::print("Allocation of ja failed.");
                error -100002;
            }
            std::copy(ja0.data(),ja0.data()+ja0.length(),ja);
            a = (mreal*) malloc(a0.length()*sizeof(mreal));
            if(!a)
            {
                mma::print("Allocation of a failed.");
                error -100003;
            }
            std::copy(a0.data(),a0.data()+a0.length(),a);
            perm = (mint*) malloc(n*sizeof(mint));
            if(!perm)
            {
                mma::print("Allocation of perm failed.");
                error -100004;
            }
            for ( i = 0; i < n; i++ )
            {
                perm[i] = i+1;
            }

            maxfct = 1;                 /* Maximum number of numerical factorizations. */
            mnum = 1;                   /* Which factorization to use. */
            msglvl = 0;                 /* Print statistical information in file */
            error = 0;                  /* Initialize error flag */

            iparm[0] = 1;               /* No solver default */
            iparm[1] = 3;               /* parallel version of nested disection */
            iparm[3] = 0;               /* No iterative-direct algorithm */
            iparm[4] = 2;               /* Write fill-in reducing permutation to perm */
            iparm[5] = 0;               /* Write solution into x */
            if(mtype==11)
            {
                iparm[9] = 13;          /* Perturb the pivot elements with 1E-iparm[9] */
            }
            else
            {
                iparm[9] = 8;           /* Perturb the pivot elements with 1E-iparm[9] */
            }
            if( (mtype!=2) || (mtype!=-2) )
            {
                iparm[10] = 0;          /* Disable scaling. Because it is slow.*/
                iparm[12] = 0;          /* Disable matching. Because it is slow.*/
            }
            else
            {
                iparm[10] = 1;          /* Enable scaling. Default for nonsymmetric matrices. Good for indefinite symmetric matrices */
                iparm[12] = 1;          /* Enable matching. Default for nonsymmetric matrices. Good for indefinite symmetric matrices */
            }
            iparm[17] = -1;             /* Report number of nonzeros in the factor LU */
            iparm[18] = -1;             /* Report Mflops for LU factorization */

            iparm[20] = 1;              /* Bunch-Kaufman pivoting */
            iparm[34] = 0;              /* 1-based indexing */

            if(error==0)
            {
                //mma::print("Initialization done.");
                initialized = 1;
            }
            else
            {
                mma::print("Initialization failed.");
                initialized = 0;
            }
        }
	}
	mint SetNonzeroValues(mma::RealTensorRef a0)
	{
		if(a)
		{
			std::copy(a0.data(),a0.data()+nnz,a);
            numfactorized = 0;
            return 0;
		}
        else
        {
            mma::print("SetNonzeroValues failed. No memory allocated.");
            numfactorized = 0;
            return -10001;
        }

	}
	mint SetIntegerParameters(mma::IntTensorRef iparm0)
	{
		std::copy(iparm0.data(),iparm0.data()+64,iparm);
		symfactorized = 0;
		numfactorized = 0;
		return 0;
	}
	mint SetTimeStamp(mreal t)
	{
		timestamp = t;
		return 0;
	}
    mint SetCheckMatrixQ(mint c)
    {
        checkmatrix = c;
        return 0;
    }
	mint SetPermutations(mma::IntTensorRef perm0)
	{
        if(perm)
        {
            std::copy(perm0.data(),perm0.data()+n,perm);
            return 0;
        }
        else
        {
            mma::print("SetPermutations failed. No memory allocated.");
            return -10000;
        }
	}
	mint FactorizeSymbolically()
	{
		if(!initialized)
		{
			mma::print("Pardiso has to be initialized before symbolic factorization.");
			error = -111;
			return error;
		}
        //mma::print("Starting symbolic factorization.");
		phase = 11;
		pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
		if(error!=0)
		{
			symfactorized = 0;
            mma::print("Pardiso reported an error in symbolic factorization.");
		}
        else
        {
            //mma::print("Symbolic factorization completed.");
            symfactorized = 1;
        }
        return error;
	}
	mint FactorizeNumerically()
	{
		if(!symfactorized)
		{
			FactorizeSymbolically();
		}
		phase = 22;
        //mma::print("Starting numeric factorization.");
		pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
		if(error!=0)
		{
            numfactorized = 0;
            mma::print("Pardiso reported an error in numeric factorization.");
		}
        else
        {
            //mma::print("Numeric factorization completed.");
            numfactorized = 1;
        }
		return error;
	}
	mma::RealTensorRef LinearSolve(mma::RealTensorRef b, mint nrhs0, mint mode)
	{
		if(!numfactorized)
		{
			FactorizeNumerically();
		}
		auto x = mma::makeVector<mreal>(b.length());
		phase = 33;
		iparm[11] = mode;
		nrhs = nrhs0;
		pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, b.data(), x.data(), &error);
		if(error!=0)
		{
			mma::print("Pardiso reported an error in solving phase.");
		}
		nrhs = 1;
		iparm[11] = 0;
		return x;
	}
	
	protected:
	;
};
;
