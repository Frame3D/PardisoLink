#include <algorithm>
#include "LTemplate.h"

extern "C" void pardiso (void *, mint *, mint *, mint *, mint *, mint *, mreal *, mint *, mint *, mint *, mint *, mint *, mint *, mreal *, mreal *, mint *);
class Pardiso
{
    mint n =0;                      /* Length of matrix */
    mint nnz =0;                    /* Number of nonzeroes */
    mint nrhs =1;                   /* Number of right hand sides */
    mint mtype =0;                  /* Matrix type */
    mint mnum = 1;                  /* Which factorization to use */
    mint error = 0;                 /* Error flag */
    mint msglvl = 0;                /* Do not pint statistical information to file */
    mint maxfct = 1;                /* Maximum number of numerical factorizations */
    mint phase = 0;                 /* Phase of pardiso */
    mint *perm = NULL;              /* Permutation */
    mint *ia = NULL;                /* Row pointers */
    mint *ja = NULL;                /* Column indices */
	mreal* a = NULL;                /* Array of nonzero values */
	mint iparm[64];                 /* Integer parameter array for controlling pardiso */

	void* pt[64];                   /* Pointer used internally by pardiso to store its data */
	mreal ddum;                     /* double dummy pointer */
	mint idum;                      /* integer dummy pointer */
	mreal timestamp = 0;            /* time stamp of the numerica factorization */

    /* a couple of flags */
	mint initialized = 0;
	mint symfactorized = 0;
	mint numfactorized = 0;
    mint checkmatrix = 1;
	
	public:
	;
	Pardiso()
	{
        mint i;
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
        else
        {
            return mma::makeVector<mint>(0,ia);
        }
	}
	mma::IntTensorRef ColumnIndices()
	{
        if(ja)
        {
            return mma::makeVector<mint>(nnz,ja);
        }
        else
        {
            return mma::makeVector<mint>(0,ja);
        }
	}
	mma::RealTensorRef NonzeroValues()
	{
        if(a)
        {
            return mma::makeVector<mreal>(nnz,a);
        }
        else
        {
            return mma::makeVector<mreal>(0,a);
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
        else
        {
            return mma::makeVector<mint>(0,perm);
        }
	}
	void Init(mma::IntTensorRef ia0, mma::IntTensorRef ja0, mma::RealTensorRef a0, mint mtype0)
	{
        if(!initialized)
        {
            mint i;
            mtype = mtype0;
            nrhs = 1;
            timestamp = 0.;
            n = ia0.length()-1;
            nnz = ia0[ia0.length()-1]-1;
            ia = (mint*) malloc(ia0.length()*sizeof(mint));
            mint anyerror = 0;
            if(!ia)
            {
                mma::print("Allocation of ia failed.");
                anyerror = 1;
            }
            std::copy(ia0.data(),ia0.data()+ia0.length(),ia);
            ja = (mint*) malloc(ja0.length()*sizeof(mint));
            if(!ja)
            {
                mma::print("Allocation of ja failed.");
                anyerror = 1;
            }
            std::copy(ja0.data(),ja0.data()+ja0.length(),ja);
            a = (mreal*) malloc(a0.length()*sizeof(mreal));
            if(!a)
            {
                mma::print("Allocation of a failed.");
                anyerror = 1;
            }
            std::copy(a0.data(),a0.data()+a0.length(),a);
            perm = (mint*) malloc(n*sizeof(mint));
            if(!perm)
            {
                mma::print("Allocation of perm failed.");
                anyerror = 1;
            }
            for ( i = 0; i < n; i++ )
            {
                perm[i] = i+1;
            }

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
            iparm[18] = 0;              /* Do not compute Mflops for LU factorization (because it is not for free) */

            iparm[20] = 1;              /* Bunch-Kaufman pivoting */
            iparm[34] = 0;              /* 1-based indexing */

            if(anyerror==0)
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
	mma::RealTensorRef LinearSolve(mma::RealTensorRef b, mint mode)
	{
		if(!numfactorized)
		{
			FactorizeNumerically();
		}
		auto x = mma::makeVector<mreal>(b.length());
		phase = 33;
		iparm[11] = mode;
		nrhs = 1;
		pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, b.data(), x.data(), &error);
		if(error!=0)
		{
			mma::print("Pardiso reported an error in solving phase.");
		}
		nrhs = 1;
		iparm[11] = 0;
		return x;
	}

    mma::RealTensorRef LinearSolveMatrix(mma::RealMatrixRef B, mint mode)
    {
        if(!numfactorized)
        {
            FactorizeNumerically();
        }
        auto X = mma::makeMatrix<mreal>(B.rows(),B.cols());
        phase = 33;
        iparm[11] = mode;
        nrhs = B.rows();
        pardiso (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, B.data(), X.data(), &error);
        if(error!=0)
        {
            mma::print("Pardiso reported an error in solving phase.");
        }
        nrhs = 1;
        iparm[11] = 0;
        return X;
    }
	
	protected:
	;
};
;
