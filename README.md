# sfh_parameter_optimizer

Generic C header-only arbitrary function parameter optimizer

If you have an arbitrary function and you need to find a series of inputs to minimize the error for the function, you can use this tool.

By default, it can somewhat work its way out of some local minimums, it is not a fully simulated annealing system, but could be modified to be such.

You just need to provide the following, per-variable
```c
	// For any time the optimizer wants to adjust some variable.
	// This needs to updated the variable, and, cascade any error calculations.
	void AdjustFunction( void * vOpaque, int iOpaque, double v );
```

And the following for the whole context:
```c
	// Provide the optimizer the system's current error.
	double ComputeError( void * vOpaque, int iOpaque );

	// Should the optimizer continue or terminate?
	int ContinueGoing( struct optimum_ctx_t * ctx );
```

You can setup a system like the following:

```c

double my_var_1, my_var_2, my_var_3

void AdjustFunction( void * vOpaque, int iOpaque, double v );
{
	*((double*)vOpaque) = v;
}

double ComputeError( void * vOpaque, int iOpaque )
{
	return my_system_error;
}

int ContinueGoing( struct optimum_ctx_t * ctx )
{
	// Run for 100 iterations
	return ctx->iteration < 100;
}

main()
{
	optimum_variable variables[4];

	variables[0] = (optimum_variable){ AdjustGlobal, &my_var_1, 0, 0.04, 1.0, 0.001 };
	variables[1] = (optimum_variable){ AdjustGlobal, &my_var_2, 0, 0.04, 1.0, 0.001 };
	variables[2] = (optimum_variable){ AdjustGlobal, &my_var_3, 0, 0.01, 0.0, 0.001 };
	variables[3] = (optimum_variable){ AdjustGlobal, &my_var_4, 0, 0.01, 0.0, 0.001 };

	optimum_ctx ctx = {
		.fnCompute = GetGlobalError,
		.fnContinueGoing = ContinueGoing,
		.vOpaque = 0,
		.iOpaque = 0,
		.nrVariables = 4,
		.variables = variables,
	};
	double er = SFHGenericOptimumOptimize( &ctx );
	return er;
}
```


