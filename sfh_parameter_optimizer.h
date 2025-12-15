#ifndef _SFH_GENERIC_OPTIMUM_H
#define _SFH_GENERIC_OPTIMIM_H

// This isn't great at finding nonlocal minimums, but it does kinda work.
// see grow and shrink constants.

// XXX TODO: We could be more celver about not re-evaluating things
// that don't need to change.

/*
	2025 CNLohr,

	This is free and unencumbered software released into the public domain.

	Anyone is free to copy, modify, publish, use, compile, sell, or
	distribute this software, either in source code form or as a compiled
	binary, for any purpose, commercial or non-commercial, and by any
	means.

	In jurisdictions that recognize copyright laws, the author or authors
	of this software dedicate any and all copyright interest in the
	software to the public domain. We make this dedication for the benefit
	of the public at large and to the detriment of our heirs and
	successors. We intend this dedication to be an overt act of
	relinquishment in perpetuity of all present and future rights to this
	software under copyright law.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
	EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
	MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
	IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
	OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
	ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
	OTHER DEALINGS IN THE SOFTWARE.

	For more information, please refer to <https://unlicense.org/>
*/

struct  optimum_ctx_t;

// Don't make this too small, otherwise you'll do poor local minimum handling.
#define GENERIC_OPT_SHRINK 0.97

// Set this to be >2 to make it easier to get out of poor local minimums.
#define GENERIC_OPT_GROW 1.9


typedef void (*GEOAdjustFunction)( void * vOpaque, int iOpaque, double v );
typedef double (*GEOComputeError)( void * vOpaque, int iOpaque );
typedef int (*GEOContinueGoing)( struct optimum_ctx_t * ctx );

typedef struct
{
	// These must be initialized by the user.
	GEOAdjustFunction fnAdjust;
	void * vOpaque;
	int iOpaque;
	double maxStep;
	double currentValue;
	double stepSize;

	// These are internal.
} optimum_variable;

struct optimum_ctx_t
{
	GEOComputeError fnCompute;
	GEOContinueGoing fnContinueGoing;
	void *vOpaque;
	int iOpaque;

	int nrVariables;
	optimum_variable * variables;

	int iteration;
	double currentError;
	double previousError;
};

typedef struct optimum_ctx_t optimum_ctx;

static inline void SFHGenericOptimumTune( optimum_ctx * ctx, optimum_variable * var )
{
	double startError = ctx->currentError;
	double startValue = var->currentValue;
	double stepSize = var->stepSize;

	var->currentValue = startValue + stepSize;
	var->fnAdjust( var->vOpaque, var->iOpaque, var->currentValue );
	double addError = ctx->fnCompute( ctx->vOpaque, ctx->iOpaque );

	var->currentValue = startValue - stepSize;
	var->fnAdjust( var->vOpaque, var->iOpaque, var->currentValue );
	double subError = ctx->fnCompute( ctx->vOpaque, ctx->iOpaque );

	// Based on start, add and sub, figure out what to do.
	if( addError > startError && subError > startError )
	{
		// Do not adjust currentValue
		var->stepSize = stepSize * GENERIC_OPT_SHRINK;
		var->currentValue = startValue;
		var->fnAdjust( var->vOpaque, var->iOpaque, startValue );
		ctx->currentError = startError;
	}
	else
	{
		// Going up is better.
		if( addError < subError )
		{
			var->currentValue = startValue + stepSize;
			var->fnAdjust( var->vOpaque, var->iOpaque, var->currentValue );
			ctx->currentError = addError;
		}
		else
		{
			ctx->currentError = subError;

			// Sub error is second, no need to adjust it back.
		}

		stepSize = stepSize * GENERIC_OPT_GROW;
		if( stepSize >= var->maxStep )
			stepSize = var->maxStep;
		var->stepSize = stepSize;
	}

}

static inline double SFHGenericOptimumOptimize( optimum_ctx * ctx )
{
	int v;
	for( v = 0; v < ctx->nrVariables; v++ )
	{
		optimum_variable * o = ctx->variables + v;
		o->fnAdjust( o->vOpaque, o->iOpaque, o->currentValue );
	}
	ctx->previousError = 1e30;
	ctx->currentError = ctx->fnCompute( ctx->vOpaque, ctx->iOpaque );
	ctx->iteration = 0;
	while( ctx->fnContinueGoing( ctx ) )
	{
		double prevError = ctx->currentError;
		ctx->iteration++;
		for( v = 0; v < ctx->nrVariables; v++ )
		{
			optimum_variable * o = ctx->variables + v;
			SFHGenericOptimumTune( ctx, o );
		}
		ctx->currentError = ctx->fnCompute( ctx->vOpaque, ctx->iOpaque );
		ctx->previousError = prevError;
	}
	return ctx->currentError;
}

#endif


