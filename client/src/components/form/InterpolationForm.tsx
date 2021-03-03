import * as React from "react";
import { ICollection, IInterpolatedParticle } from "../../interfaces";
import { interpolateStructures } from "../../utils/api";
import { FormWrapper } from "./FormWrapper";
import { StructureImage } from "../StructureImage";

export const InterpolationForm = ({
    open,
  setCollection,
  structures,
  setStructures,
  loading,
  setLoading,
}: {
  open?: boolean;
  setCollection(collection: ICollection<IInterpolatedParticle>): void;
  structures?: string[];
  setStructures?(structures: string[]): void;
  loading: boolean;
  setLoading(loading: boolean): void;
}) => {
  const [name, setName] = React.useState<string>("Interpolated");

  return (
    <FormWrapper
      open={open}
      title="Interpolate between structures"
      loading={loading}
      setLoading={setLoading}
      onSubmit={async () => {
        let data: IInterpolatedParticle[] = [];
        const validStructures = structures?.filter(Boolean);
        if (validStructures && validStructures.length > 0) {
          setLoading(true);
          data = await interpolateStructures(validStructures, 30);
        }
        setCollection({
          data,
          name,
        });
        setLoading(false);
      }}
    >
          {/* We rely on interpolation collection to be called "Interpolated" 
           <div className="form-group">
            <label htmlFor="interpolationNameInput">Name</label>
            <input
              type="text"
              className="form-control form-control-sm"
              id="interpolationNameInput"
              required
              value={name}
              onChange={(e) => setName(e.currentTarget.value)}
            />
          </div> */}
          <div className="form-group">
            <label htmlFor="inpolationStructures">Structures (newline separated)</label>
            <textarea
              className="form-control form-control-sm"
              id="inpolationStructures"
              rows={3}
              onChange={(e) => setStructures?.(e.currentTarget.value.split("\n"))}
              value={structures?.join("\n") || ""}
            />
          </div>
          <div className="text-right">
            <button className="btn btn-primary" type="submit" disabled={loading}>
              {loading ? (
                <>
                  <span className="spinner-grow spinner-grow-sm" role="status" aria-hidden="true"></span> Loading...
                </>
              ) : (
                <>Compute Interpolation</>
              )}
            </button>
          </div>
          <div className="text-center">
            {structures && structures.length > 0
              ? structures
                  .map<React.ReactNode>((structure) => (
                    <StructureImage key={structure} structure={structure} width="70px" height="70px" />
                  ))
                  .reduce((prev, curr, i) => [
                    prev,
                    <i key={i} className="fas  fa-fw fa-long-arrow-alt-right ml-2 mr-2" />,
                    curr,
                  ])
              : null}
          </div>
    </FormWrapper>
  );
};
