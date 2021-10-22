import * as React from "react";
import { ICollection, IParticleSelection } from "../../interfaces";
import { computeProjectionsWithModels } from "../../utils/api";
import { ButtonWithUpload } from "../ButtonWithUpload";
import { FormWrapper } from "./FormWrapper";

export const RecomputeEmbeddingsForm = ({
  collections,
  setCollections,
  selection,
  loading,
  setLoading,
}: {
  collections: ICollection[];
  setCollections(collections: ICollection[]): void;
  selection: IParticleSelection;
  loading: boolean;
  setLoading(loading: boolean): void;
}) => {
  const [from, setFrom] = React.useState<ICollection | null>(null);
  const [projection, setProjection] = React.useState<string | null>("all");
  const [to, setTo] = React.useState<ICollection | null>(null);

  const recomputeProjections = async () => {
    if (from?.projections && to && projection) {
      const recomputedParticles = await computeProjectionsWithModels(
        to.data,
        projection === "all" ? Object.entries(from.projections).reduce((acc, [key, value]) => ({...acc, [key]: value.model}), {}) : { [projection]: from.projections[projection].model }
      );

      const toName = (projection: string): string => `${from.name}_${projection}_base`;

      const newTo = {
        ...to,
        data: to.data.map((d, i) => ({
          ...d,
          projection: {
            ...d.projection,
            ...Object.entries(recomputedParticles.projection).reduce(
              (acc, [projection, data]) => ({
                ...acc,
                [toName(projection)]: data[i],
              }),
              {}
            ),
          },
        })),
        projections: {
          ...(to.projections || {}),
          ...Object.entries(recomputedParticles.projections || {}).reduce((acc, [projection, data]) => ({
            ...acc,
            [toName(projection)]: data
          }), {})
        }
      };

      const newFrom = {
        ...from,
        data: from.data.map((d) => ({
          ...d,
          projection: {
            ...d.projection,
            ...Object.keys(recomputedParticles.projection).reduce(
              (acc, projection) => ({
                ...acc,
                [toName(projection)]: d.projection[projection],
              }),
              {}
            ),
          },
        })),
      };
      // throw Error('TEST');

      setCollections(collections.map((c) => (c === to ? newTo : c === from ? newFrom : c)));

      setTo(newTo);
      setFrom(newFrom);
    }
  };

  const collectionsWithProjections = collections.filter(
    (c) => c.projections && Object.entries(c.projections).some(([key, value]) => value.model)
  );

  const collectionsWithEmbeddings = collections.filter(
    (c) => c.data?.[0].embedding && Object.entries(c.data?.[0].embedding).length > 0
  );

  return (
    <FormWrapper
      title="Recompute Projections"
      loading={loading}
      setLoading={setLoading}
      onSubmit={recomputeProjections}
    >
      {collectionsWithProjections.length > 0 ? (
        <>
          <div className="mb-3">
            <label htmlFor="recomputeProjectionFrom">Reference collection</label>
            <div className="input-group input-group-sm">
              <select
                className="form-control"
                id="recomputeProjectionFrom"
                value={from?.name}
                onChange={(e) => {
                  const from = e.currentTarget.value;
                  setFrom(from ? collections.find((c) => c.name === from)! : null);
                }}
              >
                <option value="">Select...</option>
                {collectionsWithProjections.map((c) => (
                  <option key={c.name} value={c.name}>{c.name}</option>
                ))}
              </select>
            </div>
          </div>
          {from ? (
            <div className="mb-3">
              <label htmlFor="recomputeProjectionProjection">Projection</label>
              <div className="input-group input-group-sm">
                <select
                  className="form-control"
                  id="recomputeProjectionProjection"
                  value={projection || undefined}
                  onChange={(e) => {
                    setProjection(e.currentTarget.value);
                  }}
                >
                  <option value="all">All</option>
                  {Object.entries(from?.projections || {})
                    .filter(([key, value]) => value.model)
                    .map(([key, value]) => (
                      <option value={key}>{key}</option>
                    ))}
                </select>
              </div>
            </div>
          ) : null}
          <div className="mb-3">
            <label htmlFor="recomputeProjectionTo">Subset collection</label>
            <div className="input-group input-group-sm">
              <select
                className="form-control"
                id="recomputeProjectionTo"
                value={to?.name}
                onChange={(e) => {
                  const to = e.currentTarget.value;
                  setTo(to ? collections.find((c) => c.name === to)! : null);
                }}
              >
                <option value="">Select...</option>
                {collectionsWithEmbeddings.map((c) => (
                  <option value={c.name}>{c.name}</option>
                ))}
              </select>
            </div>
          </div>
          <div className="text-end">
            <ButtonWithUpload
              loading={loading}
              disabled={Boolean(!from || !to || from === to)}
              text="Compute Projections"
            />
          </div>
        </>
      ) : (
        <p>At least 2 collections must be available.</p>
      )}
    </FormWrapper>
  );
};
