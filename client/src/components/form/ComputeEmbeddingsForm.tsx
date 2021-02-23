import * as React from "react";
import { DEFAULT_COLLECTION, ICollection } from "../../interfaces";
import { embedStructures } from "../../utils/api";
import { ButtonWithUpload } from "../ButtonWithUpload";
import { FormWrapper } from "./FormWrapper";

export const ComputeEmbeddingsForm = ({
  addCollection,
  loading,
  setLoading,
}: {
  addCollection(collection: ICollection): void;
  loading: boolean;
  setLoading(loading: boolean): void;
}) => {
  const [cdddInput, setCdddInput] = React.useState<string>("");

  return (
    <FormWrapper
      title="Compute Embeddings"
      loading={loading}
      setLoading={setLoading}
      onSubmit={async () => {
        const structures: string[] = [];
        const additional: { [key: string]: string[] } = {};
        cdddInput
          .split("\n")
          .filter(Boolean)
          .forEach((i) => {
            const [smiles, ...rest] = i.split(/[;,\t]/);
            structures.push(smiles);
            rest.forEach((v, i) => {
              additional[i.toString()] = [...(additional[i.toString()] || []), v];
            });
          });

        const data = await embedStructures({
          structures,
          additional,
        });
        addCollection({
          data,
          name: DEFAULT_COLLECTION,
        });
      }}
    >
      <div className="form-group">
        <label htmlFor="cdddTextarea">Structures to convert</label>
        <textarea
          className="form-control form-control-sm"
          id="cdddTextarea"
          rows={3}
          onChange={(e) => setCdddInput(e.currentTarget.value)}
          value={cdddInput}
        />
      </div>
      <div className="text-right">
        <ButtonWithUpload
          loading={loading}
          text="Compute"
          onUploadData={(value) => {
            setCdddInput(value || "");
            return true;
          }}
          onUploadResult={(value) => {
            if (value) {
              addCollection({
                data: JSON.parse(value),
                name: DEFAULT_COLLECTION,
              });
              return true;
            }
          }}
        />
      </div>
    </FormWrapper>
  );
};
