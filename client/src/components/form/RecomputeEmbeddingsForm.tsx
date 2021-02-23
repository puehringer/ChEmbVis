import * as React from "react";
import { ICollection, IParticle } from "../../interfaces";
import { computeProjections } from "../../utils/api";
import { FormWrapper } from "./FormWrapper";

// TODO: Update as necessary
export const RecomputeEmbeddingsForm = ({
  collections,
  setCollections,
  particles,
  setParticles,
  interpolatedStructures,
  setInterpolatedStructures,
  loading,
  setLoading,
}: {
  collections: ICollection[];
  setCollections(collections: ICollection[]): void;
  particles: IParticle[];
  setParticles(particles: IParticle[]): void;
  interpolatedStructures: IParticle[];
  setInterpolatedStructures(particles: IParticle[]): void;
  loading: boolean;
  setLoading(loading: boolean): void;
}) => {
  const [includeInterpolatedInProjection, setIncludeInterpolatedInProjection] = React.useState<boolean>(true);

  const recomputeEmbedding = async () => {
    const validParticles = particles.filter((p) => p.embedding);
    const interpolatedParticles = includeInterpolatedInProjection
      ? interpolatedStructures.filter((p) => p.embedding)
      : undefined;
    const data = await computeProjections(
      validParticles.map((p) => p.structure),
      validParticles.map((p) => p.embedding!),
      interpolatedParticles?.map((p) => p.structure),
      interpolatedParticles?.map((p) => p.embedding!)
    );
    // Use the data to patch the embeddings
    setParticles(
      particles.map((p, i) => ({
        ...p,
        projection: data.projections.reduce((acc, key) => {
          return { ...acc, [key]: data.projection[key][i] };
        }, {}),
      }))
    );
    setInterpolatedStructures(
      data.additional
        ? interpolatedStructures.map((p, i) => ({
            ...p,
            projection: data.projections.reduce((acc, key) => {
              return { ...acc, [key]: data.additional![key]?.[i] };
            }, {}),
          }))
        : []
    );
  };

  return (
    <FormWrapper
      title="Compute projections based on the embeddings"
      loading={loading}
      setLoading={setLoading}
      onSubmit={recomputeEmbedding}
    >
      <div className="form-check">
        <input
          className="form-check-input"
          type="checkbox"
          checked={includeInterpolatedInProjection}
          onChange={(e) => setIncludeInterpolatedInProjection(e.currentTarget.checked)}
          disabled={interpolatedStructures.length === 0}
          id="includeInterpolatedCheckbox"
        />
        <label className="form-check-label" htmlFor="includeInterpolatedCheckbox">
          Include interpolated
        </label>
      </div>
      <div className="text-right mt-2">
        <button className="btn btn-primary" type="submit" disabled={loading || particles.length <= 1}>
          {loading ? (
            <>
              <span className="spinner-grow spinner-grow-sm" role="status" aria-hidden="true"></span> Loading...
            </>
          ) : (
            <>Recompute Embedding</>
          )}
        </button>
      </div>
    </FormWrapper>
  );
};
