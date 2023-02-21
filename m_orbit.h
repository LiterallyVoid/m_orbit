// Numerical stability status:
// - Elliptical, with zero velocity at apoapsis: stable
// - Parabolic, intersecting the origin: UNSTABLE
// - Hyperbolic, with turn angle near π (e≈1): UNSTABLE

#ifndef M_ORBIT_H
#define M_ORBIT_H

#ifndef MORB_NOSTD

#include <math.h>

#define MORB_FABS fabs

#define MORB_SQ(x)  ((x)*(x))
#define MORB_SQRT sqrt
#define MORB_POW pow

#define MORB_SIN sin
#define MORB_COS cos

#define MORB_LOG10 log10

#define MORB_TAN tan
#define MORB_ATAN atan
#define MORB_ATAN2 atan2

#define MORB_SINH sinh
#define MORB_COSH cosh

#define MORB_TANH tanh
#define MORB_ATANH atanh

#define MORB_PI M_PI

#define MORB_VEC_EL double

#endif

typedef MORB_VEC_EL morbEl;

typedef enum {
	MORB_ORBIT_ELLIPTIC,
	MORB_ORBIT_PARABOLIC,
	MORB_ORBIT_HYPERBOLIC,
} morbKind;

typedef struct {
	morbEl v[4];
} morbVec;

typedef struct {
	morbKind kind;

	morbEl energy;

	morbEl eccentricity;

	// For eccentricites very close to 1, we assume that the eccentricity is
	// instead 1-1e6 and squash the semiminor axis of the ellipse to compensate.
	morbEl semi_minor_scale;

	// Parameters for converting Equatorial <-> Perifocal coordinates.
	morbEl longitude_of_ascending_node;
	morbEl inclination;

	morbEl argument_of_periapsis;

	// Monotonic location of the body along the orbit.  For parabolic orbits,
	// mean anomaly is undefined, and is substituted for time since periapsis.
	morbEl mean_anomaly;
} morbElements;

typedef struct {
	morbVec position;
	morbVec velocity;
} morbVectors;

// `mu` is the gravitational constant times the orbiting center mass.
morbVectors morbElementsToVectorsAtTrueAnomaly(morbElements elements, morbEl mu, morbEl true_anomaly);

morbVectors morbElementsToVectors(morbElements elements, morbEl mu, morbEl tolerance);
morbElements morbVectorsToElements(morbVectors vectors, morbEl mu);

morbEl morbMeanAnomalyPerTime(morbElements elements, morbEl mu);

#endif // M_ORBIT_H

#if defined(M_ORBIT_IMPLEMENTATION)

static morbVec morbAdd(morbVec a, morbVec b) {
	return (morbVec) {
		{
			a.v[0] + b.v[0],
			a.v[1] + b.v[1],
			a.v[2] + b.v[2],
			a.v[3] + b.v[3],
		},
	};
}

static morbVec morbSub(morbVec a, morbVec b) {
	return (morbVec) {
		{
			a.v[0] - b.v[0],
			a.v[1] - b.v[1],
			a.v[2] - b.v[2],
			a.v[3] - b.v[3],
		},
	};
}

static morbVec morbMul(morbVec a, morbVec b) {
	return (morbVec) {
		{
			a.v[0] * b.v[0],
			a.v[1] * b.v[1],
			a.v[2] * b.v[2],
			a.v[3] * b.v[3],
		},
	};
}

static morbVec morbScale(morbVec vec, morbEl scale) {
	return (morbVec) {
		{
			vec.v[0] * scale,
			vec.v[1] * scale,
			vec.v[2] * scale,
			vec.v[3] * scale,
		}
	};
}

static morbVec morbSclAddSclAdd(morbVec initial, morbVec a, morbEl a_scale, morbVec b, morbEl b_scale) {
	return morbAdd(
		initial,
		morbAdd(
			morbScale(a, a_scale),
			morbScale(b, b_scale)
		)
	);
}

static morbEl morbLength3(morbVec vec) {
	return MORB_SQRT(
		vec.v[0] * vec.v[0] +
		vec.v[1] * vec.v[1] +
		vec.v[2] * vec.v[2]
	);
}

static morbVec morbNormalize3(morbVec vec) {
	morbEl length = morbLength3(vec);

	if (length == 0.0) {
		return (morbVec) {
			{
				0.0,
				0.0,
				1.0,
				0.0,
			},
		};
	}

	return morbScale(vec, 1.0 / length);
}

static morbEl morbDot(morbVec a, morbVec b) {
	return
		a.v[0] * b.v[0] +
		a.v[1] * b.v[1] +
		a.v[2] * b.v[2] +
		a.v[3] * b.v[3];
}

// Equivalent to morbCross3({0, 0, 1}, vec)
static morbVec morbCross2(morbVec vec) {
	return (morbVec) {
		{
			-vec.v[1],
			vec.v[0],
			0.0,
			0.0,
		},
	};
}

static morbVec morbCross3(morbVec a, morbVec b) {
	return (morbVec) {
		{
			a.v[1] * b.v[2] - a.v[2] * b.v[1],
			a.v[2] * b.v[0] - a.v[0] * b.v[2],
			a.v[0] * b.v[1] - a.v[1] * b.v[0],
			0.0,
		},
	};
}

static morbEl morbTrueAnomalyFromMeanAnomalyElliptic(morbEl eccentricity, morbEl mean_anomaly, morbEl tolerance) {
	const morbEl e = eccentricity;

	const morbEl M = mean_anomaly;
	morbEl E = mean_anomaly - e * MORB_SIN(mean_anomaly);

	for (int i = 0; i < 30; i++) {
		const morbEl D_M = M - (E - e * MORB_SIN(E));
		const morbEl D_E = D_M / (1 - (e * MORB_COS(E)));

		E += D_E;

		if (MORB_FABS(D_E) < tolerance) {
			break;
		}
	}

	const morbEl beta = e / (1 + MORB_SQRT(1 - (e * e)));
	const morbEl true_anomaly = E + 2 * MORB_ATAN(
		(beta * MORB_SIN(E)) / (1 - beta * MORB_COS(E))
	);

	return true_anomaly;
}

static morbEl morbTrueAnomalyFromMeanAnomalyHyperbolic(morbEl eccentricity, morbEl mean_anomaly, morbEl tolerance) {
	const morbEl e = eccentricity;

	const morbEl M = mean_anomaly;
	morbEl F = 0.1;
	if (M > 20) {
		F = MORB_LOG10(M) * 0.7 + 3;
	}

	for (int i = 0; i < 30; i++) {
		const morbEl D_M = e * MORB_SINH(F) - F - M;
		const morbEl D_F = -D_M / (e * MORB_COSH(F) - 1.0);

		F += D_F;

		if (MORB_FABS(D_F) < tolerance) {
			break;
		}
	}

	morbEl true_anomaly = 2.0 * MORB_ATAN(
		MORB_TANH(F / 2.0) / MORB_SQRT((e - 1) / (e + 1))
	);

	return true_anomaly;
}

morbEl morbMeanAnomalyPerTime(morbElements elements, morbEl mu) {
	if (elements.kind == MORB_ORBIT_PARABOLIC) {
		return MORB_SQ(mu) / MORB_POW(elements.energy, 3.0);
	} else if (elements.kind == MORB_ORBIT_ELLIPTIC) {
		return (MORB_SQ(mu) / MORB_POW(elements.energy, 3.0)) * MORB_POW(1 - MORB_SQ(elements.eccentricity), 3.0 / 2.0);
	} else {
		return MORB_SQ(mu) / MORB_POW(elements.energy, 3.0) * MORB_POW(MORB_SQ(elements.eccentricity) - 1.0, 3.0 / 2.0);
	}
}

morbVectors morbElementsToVectorsAtTrueAnomaly(morbElements elements, morbEl mu, morbEl true_anomaly) {
	// Transform into the perifocal frame by:
	// - 1) rotate by the argument of periapsis around the Z axis
	// - 2) rotate by the inclination around the ascending node
	// - 3) rotate by the longitude of ascending node (Ω) around the Z axis

	const morbEl sin_O = MORB_SIN(elements.longitude_of_ascending_node);
	const morbEl cos_O = MORB_COS(elements.longitude_of_ascending_node);

	const morbEl sin_i = MORB_SIN(elements.inclination);
	const morbEl cos_i = MORB_COS(elements.inclination);

	const morbEl sin_w = MORB_SIN(elements.argument_of_periapsis);
	const morbEl cos_w = MORB_COS(elements.argument_of_periapsis);

	morbVec perifocal_x = {
		{
			cos_O * cos_w - sin_O * cos_i * sin_w,
			sin_O * cos_w + cos_O * cos_i * sin_w,
			sin_i * sin_w,
			0.0,
		},
	};

	morbVec perifocal_y = {
		{
			-cos_O * sin_w - sin_O * cos_i * cos_w,
			-sin_O * sin_w + cos_O * cos_i * cos_w,
			sin_i * cos_w,
			0.0,
		},
	};

	morbVec position, velocity;

	if (elements.kind == MORB_ORBIT_HYPERBOLIC) {
		morbEl semi_major = MORB_SQ(elements.energy) / mu / (MORB_SQ(elements.eccentricity) - 1.0);

		morbEl p = MORB_SQ(elements.energy) / mu;
		morbEl r = semi_major * (MORB_SQ(elements.eccentricity) - 1.0) / (1 + elements.eccentricity * MORB_COS(true_anomaly));

		position = (morbVec) {
			{
				MORB_COS(true_anomaly) * r,
				MORB_SIN(true_anomaly) * r,
			},
		};

		morbEl mu_over_h = MORB_SQRT(mu / p);

		velocity = (morbVec) {
			{
				mu_over_h * -MORB_SIN(true_anomaly),
				mu_over_h * (elements.eccentricity + MORB_COS(true_anomaly)),
			},
		};
	} else {
		morbEl p = MORB_SQ(elements.energy) / mu;
		morbEl r = p / (1.0 + elements.eccentricity * MORB_COS(true_anomaly));

		position = (morbVec) {
			{
				MORB_COS(true_anomaly) * r,
				MORB_SIN(true_anomaly) * r * elements.semi_minor_scale,
			},
		};

		morbEl mu_over_h = MORB_SQRT(mu / p);

		velocity = (morbVec) {
			{
				mu_over_h * -MORB_SIN(true_anomaly),
				mu_over_h * (elements.eccentricity + MORB_COS(true_anomaly)) * elements.semi_minor_scale,
			},
		};
	}

	return (morbVectors) {
		morbSclAddSclAdd((morbVec) {{ 0, 0, 0, 1 }}, perifocal_x, position.v[0], perifocal_y, position.v[1]),
		morbSclAddSclAdd((morbVec) {{ 0, 0, 0, 0 }}, perifocal_x, velocity.v[0], perifocal_y, velocity.v[1]),
	};
}

morbEl morbMeanAnomalyToTrueAnomaly(morbElements elements, morbEl mean_anomaly, morbEl tolerance) {
	if (elements.kind == MORB_ORBIT_ELLIPTIC) {
		return morbTrueAnomalyFromMeanAnomalyElliptic(elements.eccentricity, mean_anomaly, tolerance);
	} else if (elements.kind == MORB_ORBIT_PARABOLIC) {
		morbEl z = MORB_POW(3 * mean_anomaly + MORB_SQRT(1.0 + MORB_SQ(3 * mean_anomaly)), 1.0 / 3.0);
		return MORB_ATAN(z - 1.0 / z) * 2.0;
	} else {
		return morbTrueAnomalyFromMeanAnomalyHyperbolic(elements.eccentricity, mean_anomaly, tolerance);
	}
}

morbVectors morbElementsToVectors(morbElements elements, morbEl mu, morbEl tolerance) {
	return morbElementsToVectorsAtTrueAnomaly(elements, mu, morbMeanAnomalyToTrueAnomaly(elements, elements.mean_anomaly, tolerance));
}

morbElements morbVectorsToElements(morbVectors vectors, morbEl mu) {
	morbElements elements;

	// This points "up", with a magnitude depending on the energy of this orbit.
	const morbVec specific_angular_momentum = morbCross3(vectors.position, vectors.velocity);

	const morbVec specific_angular_momentum_normalized = morbNormalize3(specific_angular_momentum);

	elements.energy = morbLength3(specific_angular_momentum);

	//const morbEl semi_major = -(mu * morbLength3(vectors.position)) / (morbLength3(vectors.position) * MORB_POW(morbLength3(vectors.velocity), 2.0) - 2 * mu);

	// This points at the ascending node relative to the equatorial ecliptic, on the +Z axis.
	morbVec ascending_node_vector = morbCross2(specific_angular_momentum);

	// This points at the periapsis.
	const morbVec eccentricity_vector = morbSub(
		morbScale(morbCross3(vectors.velocity, specific_angular_momentum), 1.0 / mu),
		morbNormalize3(vectors.position)
	);

	//||r ^ v||^2 / (1 - ||v ^ (r ^ v) / μ - normalize(r)||^2)

	elements.eccentricity = morbLength3(eccentricity_vector);

	// The angle of the ascending node, around the vertical axis.
	elements.longitude_of_ascending_node = MORB_ATAN2(ascending_node_vector.v[1], ascending_node_vector.v[0]);

	// Re-normalize ascending_node_vector, making sure to preserve the angle
	// already calculated for the ascending node longitude.
	ascending_node_vector = (morbVec) {
		{
			MORB_COS(elements.longitude_of_ascending_node),
			MORB_SIN(elements.longitude_of_ascending_node),
			0.0,
			0.0,
		},
	};

	// This points to the top of the orbit, projected down onto the equatorial
	// ecliptic (i.e. this vector always has no Z component.) We'll need the
	// inclined variant to calculate the argument of periapsis later.
	const morbVec ascent_vector_equatorial = morbCross2(ascending_node_vector);

	// The inclination is how much the orbit is tilted around the ascending node
	// --- between 0 and π/2 for a counter-clockwise orbit, or between π/2 and π
	// for a clockwise orbit.
	elements.inclination = MORB_ATAN2(
		morbDot(specific_angular_momentum, ascending_node_vector),
		morbDot(specific_angular_momentum_normalized, ascent_vector_equatorial)
	);

	const morbVec ascent_vector_inclined = morbNormalize3(morbCross3(specific_angular_momentum_normalized, ascending_node_vector));

	// A ω (agument of periapsis) of zero means that the periapsis is aligned
	// with the ascending node. Increasing ω rotates the periapsis
	// counter-clockwise relative to the inclination plane.
	elements.argument_of_periapsis = MORB_ATAN2(
		morbDot(eccentricity_vector, ascent_vector_inclined),
		morbDot(eccentricity_vector, ascending_node_vector)
	);

	// True anomaly around the center.
	const morbEl true_anomaly = MORB_ATAN2(
		morbDot(vectors.position, ascent_vector_inclined),
		morbDot(vectors.position, ascending_node_vector)
	) - elements.argument_of_periapsis;

	elements.semi_minor_scale = 1.0;

	if (elements.eccentricity < (1.0 - 1e-6)) {
		elements.kind = MORB_ORBIT_ELLIPTIC;
	} else if (elements.eccentricity > (1.0 + 1e-9)) {
		elements.kind = MORB_ORBIT_HYPERBOLIC;
	} else {
		const morbEl E = (MORB_SQ(morbLength3(vectors.velocity)) * 0.5 - mu / morbLength3(vectors.position));
		const morbEl a = mu / (-2.0 * E);
		const morbEl p = a * (1 - MORB_SQ(elements.eccentricity));

		if (a >= 0.0 && a < 1e+30 && p <= 1e-6) {
			elements.kind = MORB_ORBIT_ELLIPTIC;

			const morbEl b_orig = a * MORB_SQRT(1 - MORB_SQ(elements.eccentricity));
			elements.eccentricity = 1.0 - 1e-6;
			const morbEl b_new = a * MORB_SQRT(1 - MORB_SQ(elements.eccentricity));

			elements.semi_minor_scale = b_orig / b_new;

			elements.energy = MORB_SQRT(a * (1 - MORB_SQ(elements.eccentricity)));
		} else {
			elements.kind = MORB_ORBIT_PARABOLIC;
		}
	}

	// ...and thus begins the part that I don't understand and instead copied
	// from https://orbital-mechanics.space/. Seriously, check that out! It's
	// awesome!
	if (elements.kind == MORB_ORBIT_ELLIPTIC) {
		const morbEl eccentric_anomaly = MORB_ATAN(
			MORB_SQRT((1 - elements.eccentricity) / (1 + elements.eccentricity))
			* MORB_TAN(true_anomaly / 2.0)
		) * 2.0;

		const morbEl mean_anomaly = eccentric_anomaly -
			(elements.eccentricity * MORB_SQRT(1 - (elements.eccentricity * elements.eccentricity)) * MORB_SIN(true_anomaly))
			/ (1.0 + elements.eccentricity * MORB_COS(true_anomaly));

		elements.mean_anomaly = mean_anomaly;
	} else if (elements.kind == MORB_ORBIT_PARABOLIC) {
		const morbEl t_term1 = MORB_TAN(true_anomaly / 2.0) / 2.0;
		const morbEl t_term2 = MORB_POW(MORB_TAN(true_anomaly / 2.0), 3.0) / 6.0;

		elements.mean_anomaly = t_term1 + t_term2;
	} else {
		const morbEl f = MORB_ATANH(
			MORB_SQRT((elements.eccentricity - 1.0) / (elements.eccentricity + 1.0)) * MORB_TAN(true_anomaly / 2.0)
		) * 2.0;

		elements.mean_anomaly = elements.eccentricity * MORB_SINH(f) - f;
	}

	return elements;
}

#endif // M_ORBIT_IMPLEMENTATION

#ifdef M_ORBIT_TEST

#include <stdio.h>

void printVectors(morbVectors vectors) {
	printf("\tposition: (%f,%f,%f)\tvelocity: (%f,%f,%f)\n",
		vectors.position.v[0],
		vectors.position.v[1],
		vectors.position.v[2],

		vectors.velocity.v[0],
		vectors.velocity.v[1],
		vectors.velocity.v[2]
	);
}

void printElements(morbElements elements) {
	printf("\th: %f\te: %f %s\tΩ: %f\n\ti: %f\tω: %f\tM: %f\n",
		elements.energy,
		elements.eccentricity,
		elements.kind == MORB_ORBIT_PARABOLIC ? "P" :
		elements.kind == MORB_ORBIT_HYPERBOLIC ? "H" :
		"E",
		elements.longitude_of_ascending_node,
		elements.inclination,
		elements.argument_of_periapsis,
		elements.mean_anomaly
	);
}

int main(int argc, char **argv) {
	const morbEl mu = 1.0;

	morbVectors cases[] = {
		{
			{ 0, 0, 0, 1 },
			{ 0, 0, 0, 0 },
		},
		{
			{ 0.0, 1.0, 0, 1 },
			{ 1e-7, 0.0, 0, 0 },
		},
		{
			{ 0.0, -1.0, 0, 1 },
			{ -1.0, 0.0, 0, 0 },
		},
		{
			{ 0.5, 0, 0, 1 },
			{ 0, 2.0, 0, 0 },
		},
		{
			{ 1, 0, 0, 1 },
			{ 0, 0.3, 0, 0 },
		},
		{
			{ 1, 0, 0, 1 },
			{ 0, 1, 0, 0 },
		},
		{
			{ 1, 0, 0, 1 },
			{ 0, MORB_SQRT(2) - 1e-9, 0, 0 },
		},
		{
			{ 1, 0, 0, 1 },
			{ 0, MORB_SQRT(2) - 1e-10, 0, 0 },
		},
		{
			{ 1, 0, 0, 1 },
			{ 0, MORB_SQRT(2), 0, 0 },
		},
		{
			{ 1, 0, 0, 1 },
			{ 0, MORB_SQRT(2) + 1e-10, 0, 0 },
		},
		{
			{ 1, 0, 0, 1 },
			{ 0, MORB_SQRT(2) + 1e-9, 0, 0 },
		},
		{
			{ 1, 0, 0, 0 },
			{ 0, 2.0, 0, 0 },
		},
		{
			{ 1, 0, 0, 0 },
			{ 0, 20.0, 0, 0 },
		},
		{
			{ 1, 0, 0, 0 },
			{ 0, 2000.0, 0, 0 },
		},
		{
			{ 1, 0, 0, 0 },
			{ 0, 20000.0, 0, 0 },
		},
		{
			{ 1, 0, 0, 0 },
			{ MORB_SQRT(2), 0, 0, 0 },
		},
		{
			{ 1, 0, 0, 0 },
			{ 5, 0, 0, 0 },
		},
		{
			{ 1, 0, 0, 0 },
			{ 5, 1e-8, 0, 0 },
		},
		{
			{ 1, 0, 0, 0 },
			{ 5, 1e-6, 0, 0 },
		},
		{
			{ 1, 0, 0, 0 },
			{ 5, 1e-4, 0, 0 },
		},
		{
			{ 1, 0, 0, 0 },
			{ 5, 1e-2, 0, 0 },
		},
		{
			{ 1, 0, 0, 0 },
			{ 5, 1e-1, 0, 0 },
		},
		{
			{ 1, 0, 0, 0 },
			{ 5, 1, 0, 0 },
		},
		{
			{ 1, 0, 0, 0 },
			{ 5, 5, 0, 0 },
		},
	};

	for (int i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
		printf("=== CASE %d\n", i);

		printVectors(cases[i]);

		morbElements elements = morbVectorsToElements(cases[i], mu);

		printf("calc els:\n");
		printElements(elements);

		printf("calc els -> roundtrip:\n");

		morbVectors vectors = morbElementsToVectors(elements, mu, 1e-14);
		printVectors(vectors);

		float anom_per_time = morbMeanAnomalyPerTime(elements, mu);

		printf("\tM/τ: %f, period (ellip.): %f\n", anom_per_time, (MORB_PI * 2.0) / anom_per_time);

		morbEl a = (MORB_SQ(elements.energy) / 1.0) / (1 - MORB_SQ(elements.eccentricity));
		printf("semimajor: %f\n", a);
		morbEl b = a * MORB_SQRT(1 - MORB_SQ(elements.eccentricity));
		printf("semiminor: %f\n", b);

		printf("at τ=0.001:\n");
		elements.mean_anomaly += anom_per_time * 0.001;
		printVectors(morbElementsToVectors(elements, mu, 1e-14));

		morbVec velocity_derived = morbScale(morbSub(morbElementsToVectors(elements, mu, 1e-14).position, vectors.position), 1.0 / 0.001);

		printf("velocity: %f,%f,%f -- difference to base: %f, difference to calculated: %f\n",
			velocity_derived.v[0],
			velocity_derived.v[1],
			velocity_derived.v[2],

			morbLength3(morbSub(velocity_derived, cases[i].velocity)),
			morbLength3(morbSub(velocity_derived, vectors.velocity)));

		printf("at τ=0.01:\n");
		elements.mean_anomaly += anom_per_time * 0.009;
		printVectors(morbElementsToVectors(elements, mu, 1e-14));

		printf("at τ=0.1:\n");
		elements.mean_anomaly += anom_per_time * 0.09;
		printVectors(morbElementsToVectors(elements, mu, 1e-14));

		printf("at τ=1.0:\n");
		elements.mean_anomaly += anom_per_time * 0.9;
		printVectors(morbElementsToVectors(elements, mu, 1e-14));

		printf("\n\n");
	}
}

#endif // M_ORBIT_TEST
