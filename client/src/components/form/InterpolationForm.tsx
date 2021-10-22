import * as React from "react";
import { ICollection, IInterpolatedParticle } from "../../interfaces";
import { interpolateStructures } from "../../utils/api";
import { FormWrapper } from "./FormWrapper";
import { StructureImage } from "../StructureImage";
import { useNameInput } from "../../utils/hooks";

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
  const [name, setName, nameInput] = useNameInput("interpolationNameInput", "Interpolated");

  return (
    <FormWrapper
      open={open}
      title="Interpolate between Structures"
      loading={loading}
      setLoading={setLoading}
      onSubmit={async () => {
        const validStructures = structures?.filter(Boolean);
        if (validStructures && validStructures.length > 0) {
          setLoading(true);
          const serverCollection = await interpolateStructures(validStructures, 100);
          setCollection({
            name,
            ...serverCollection,
          });
          setName("");
        }
        setLoading(false);
      }}
    >
          {nameInput}
          <div className="mb-3">
            <label htmlFor="inpolationStructures">Structures (newline separated)</label>
            <textarea
              className="form-control form-control-sm"
              id="inpolationStructures"
              rows={3}
              onChange={(e) => setStructures?.(e.currentTarget.value.split("\n"))}
              value={structures?.join("\n") || ""}
            />
          </div>
          <div className="text-center">
            {structures && structures.length > 0
              ? structures
                  .map<React.ReactNode>((structure) => (
                    <StructureImage key={structure} structure={structure} width="70px" height="70px" />
                  ))
                  .reduce((prev, curr, i) => [
                    prev,
                    <i key={i} className="fas  fa-fw fa-long-arrow-alt-end ms-2 me-2" />,
                    curr,
                  ])
              : null}
          </div>
          <div className="text-end">
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
    </FormWrapper>
  );
};
